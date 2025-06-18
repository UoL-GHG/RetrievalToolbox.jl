"""
Pretty printing for OCO solar model types

$(SIGNATURES)
"""
function show(io::IO, ::MIME"text/plain", sm::OCOHDFSolarModel)

    println(io, "OCO HDF solar model from file: $(sm.file_name)")

end

"""
Brief pretty printing for OCO solar model types

$(SIGNATURES)
"""
function show(io::IO, sm::OCOHDFSolarModel)

    print(io, "OCOHDFSolarModel: $(sm.file_name)")

end



"""
Calculates the down-sampled solar spectrum at the
high-resolution grid specified within the spectral
window `swin` and saves it in `rt.hires_solar`.

$(TYPEDSIGNATURES)

# Details

In this function, the solar Doppler shift is considered via
sampling the original solar spectrum at acccording wavelengths.
The Doppler factor is defined as the relative velocity between
the observation point on Earth (for Earth-Backscatter spectra),
and a negative sign indicates the observation point moving closer
to the sun, as a fraction of the speed of light in vacuum.
"""
function calculate_solar_irradiance!(
    rt::AbstractRTMethod,
    swin::AbstractSpectralWindow,
    solar_model::AbstractSolarModel;
    doppler_factor=nothing
)

    # If no Doppler factor is provided, we could
    # calculate it right here using rt.scene.location
    # and rt.scene.time.
    if isnothing(doppler_factor)

        doppler_factor = calculate_solar_doppler_shift(
            rt.scene
        )

    end

    # Doppler effect depends on the spectral unit.

    if solar_model.ww_unit isa Unitful.LengthUnits
        # Very cheeky way of moving from λ -> λ * (1 + Doppler)
        @turbo for i in eachindex(solar_model.ww)
            solar_model.ww[i] *= (1 + doppler_factor)
        end
    elseif solar_model.ww_unit isa Unitful.WavenumberUnits
        # Very cheeky way of moving from ν -> ν / (1 + Doppler)
        @turbo for i in eachindex(solar_model.ww)
            solar_model.ww[i] /= (1 + doppler_factor)
        end
    end

    #=
    We must sample the solar spectrum at Doppler-influenced
    wavelengths and produce following result:
        y = transmittance(λ * (1 + Doppler)) * continuum(λ * (1 + Doppler))
    and in a final step, multiply by the (per-wavelength) scaler
        y = y * solar_scaler
    =#

    if solar_model isa OCOHDFSolarModel
        # Sample the solar spectrum at our Doppler-influenced
        # retrieval wavelength grid
        # (and calculate continuum * transmittance in one step)
        pwl_value_1d_axb!(
            solar_model.ww,
            solar_model.transmittance,
            solar_model.continuum,
            swin.ww_grid,
            rt.hires_solar.I,
        )

    elseif solar_model isa TSISSolarModel
        # Sample the TSIS irradiance (includes both transmittance
        # and continuum) and store in hires_solar.I

        pwl_value_1d!(
            solar_model.ww,
            solar_model.irradiance,
            swin.ww_grid,
            rt.hires_solar.I
        )

    elseif solar_model isa UoLFPSolarModel
        # Sample the solar spectrum at our Doppler-influenced
        # retrieval wavelength grid
        # (and calculate continuum * transmittance in one step)
        pwl_value_1d_axb!(
            solar_model.ww,
            solar_model.transmittance,
            solar_model.continuum,
            swin.ww_grid,
            rt.hires_solar.I,
        ) 

    end

    # Apply the solar scaler
    @turbo for i in eachindex(rt.hires_solar.I)
        rt.hires_solar.S[i,1] *= rt.solar_scaler[i]
    end

    if solar_model.ww_unit isa Unitful.LengthUnits
        # Very cheeky way of moving back from λ * (1 + Doppler) -> λ
        @turbo for i in eachindex(solar_model.ww)
            solar_model.ww[i] /= (1 + doppler_factor)
        end
    elseif solar_model.ww_unit isa Unitful.WavenumberUnits
        # Very cheeky way of moving back from ν * (1 + Doppler) -> ν
        @turbo for i in eachindex(solar_model.ww)
            solar_model.ww[i] *= (1 + doppler_factor)
        end
    end



end

"""
    Converts solar model data from "ph/s/m^2/µm" to "W/m^2/µm"
"""
function convert_solar_model_to_W!(s::OCOHDFSolarModel)

    if s.irradiance_unit == u"ph/s/m^2/µm"
        @debug "Solar model in units of ph/s/m2/µm - converting to W/m2/µm!"
        @views s.continuum[:] .*= ustrip.(Ref(u"W"),
            1.0u"s^-1" .* SPEED_OF_LIGHT ./ (s.ww[:] .* u"µm") .* PLANCK
        )

        s.irradiance_unit = u"W/m^2/µm"

    else

        @warn "Units already in W/m^2/µm - skipping!"

    end

end

"""
    Converts solar model data from "W/m^2/nm" to "ph/s/m^2/µm"
"""
function convert_solar_model_to_photons!(s::TSISSolarModel)

    if s.irradiance_unit == u"W/m^2/nm"
        @debug "Solar model in units of W/m^2/nm - converting to ph/s/m^2/µm!"
        @views s.irradiance[:] .*= ustrip.(u"m^-2 * µm^-1", # We want this in 1/m2 1/µm
        s.irradiance_unit * 1.0u"s" ./ SPEED_OF_LIGHT .* (s.ww[:] .* u"nm") ./ PLANCK
        )

        s.irradiance_unit = u"ph/s/m^2/µm"

    else

        @warn "Units already in ph/s/m^2/µm - skipping!"

    end

    return true

end


"""
Reads a JPL/OCO-type solar model HDF5 file and returns
a `OCOHDFSolarModel` object.

$(TYPEDSIGNATURES)
"""
function OCOHDFSolarModel(
    filename::String,
    band_number::Integer;
    spectral_unit=:Wavelength
    )

    @assert isfile(filename) "File $(filename) is not a regular file!"

    @debug "Opening up Solar HDF file $(filename)"
    h5 = h5open(filename, "r")

    h5g_c = h5["Solar/Continuum/Continuum_$(band_number)"]
    h5g_a = h5["Solar/Absorption/Absorption_$(band_number)"]

    # Grab Absorption and Continuum
    # (values are hard-coded in the files, but not the units, so
    #  these could theoretically change in the future)

    if spectral_unit == :Wavelength

        # This is in microns
        ww_unit = u"µm"
        rad_unit = u"ph/s/m^2/µm"

        cont_ww = (1e4 ./ h5g_c["wavenumber"][:][end:-1:1])
        trans_ww = (1e4 ./ h5g_a["wavenumber"][:][end:-1:1])
        # Continuum-level values are provided in units of
        # ph/s/m2/µm, hence no further conversion between spectral radiance
        # per wavelength and spectral radiance per wavenumber is necessary.
        cont_val = h5g_c["spectrum"][:][end:-1:1]
        trans_val = h5g_a["spectrum"][:][end:-1:1]

    elseif spectral_unit == :Wavenumber

        ww_unit = u"cm^-1"
        rad_unit = u"W/m^2/cm^-1"

        # Need microns for calculations/conversions
        cont_microns = (1e4 ./ h5g_c["wavenumber"][:])

        # This is in cm^-1
        cont_ww = h5g_c["wavenumber"][:]
        trans_ww = h5g_a["wavenumber"][:]

        # Here, we need to convert from ph/s/m2/µm into
        # W/m2/cm^-1

        cont_val = h5g_c["spectrum"][:]
        # 1) convert ph/s into W
        @views cont_val[:] .*= ustrip.(Ref(u"W"),
            1.0u"s^-1" .* SPEED_OF_LIGHT ./ (cont_microns .* u"µm") .* PLANCK
        )

        #2) convert  W/m2/µm into W/m2/cm^-1
        @views cont_val[:] ./= (1e4 ./ cont_microns) .^ 2

        trans_val = h5g_a["spectrum"][:]

    end


    # Build a continuum polynomial through fitting the data
    cont_poly = Polynomials.fit(
        cont_ww,
        cont_val,
        4, #length(cont_val) - 1
    )
    # Evaluate the polynomial at all spectral window wavelengths
    # (this is reasonable since the continuum is so smooth)
    cont_resampled = cont_poly.(trans_ww)

    # Close up HDF file, all done
    close(h5)

    # Return solar model object
    return OCOHDFSolarModel(
        filename,
        band_number,
        trans_ww,
        trans_val,
        cont_resampled,
        ww_unit, # Wavelength unit
        rad_unit# Radiance unit
    )

end

"""
    Reads the full TSIS file into memory
"""
function TSISSolarModel(
    fname::String;
    spectral_unit=:Wavelength
    )

    irradiance = ncread(fname, "SSI")
    # The TSIS file has units such as
    # W m-2 nm-1, so we need to add the caret
    # such that Unitful understands them

    # Read the irradiance unit from the NetCDF file
    irradiance_u  = replace(ncgetatt(fname, "SSI", "units"), "-" => "^-", " " => "*")
    irradiance_unit = uparse(irradiance_u)

    # Read the wavelength unit from the NetCDF file
    wavelength = ncread(fname, "Vacuum Wavelength")
    wavelength_unit = uparse(ncgetatt(fname, "Vacuum Wavelength", "units"))

    # Turn wavelength into microns
    ww = ustrip.(Ref(u"µm"), wavelength * wavelength_unit)
    ww_unit = u"µm"

    # Convert wavelength to µm
    irradiance = ustrip.(Ref(u"W/m^2/µm"), irradiance * irradiance_unit)


    if spectral_unit == :Wavelength
        # Nothing to do here

    elseif spectral_unit == :Wavenumber

        # Turn irradiance into W/m^2/cm^-1 and reverse
        # array order to make them in increasing wavenumbers
        irradiance ./= (1e4 ./ ww) .^ 2
        irradiance = irradiance[end:-1:1]
        # Set the new irradiance units
        irradiance_unit = u"W/m^2/cm^-1"

        # Turn µm into cm^-1 and reverse array
        ww = 1e4 ./ ww[end:-1:1]
        ww_unit = u"cm^-1"

    else
        @error "Unknown parameter for spectral_unit: $(spectral_unit)"
    end

    # Create solar model object
    return TSISSolarModel(
        fname,
        ww,
        ww_unit,
        irradiance,
        irradiance_unit
    )

end

"""
Reads a Fraunhofer solar line list HDF5 file and returns
a `UoLFPSolarModel` object.

"""

function UoLFPSolarModel(
    filename::String,
    swin::AbstractSpectralWindow
    )

    @assert isfile(filename) "File $(filename) is not a regular file!"

    @debug "Opening up Solar HDF file $(filename)"
    h5 = h5open(filename, "r")

    #mol_mass = h5["molecular_mass"][:]
    freq = h5["freq"][:]
    stren = h5["stren"][:]
    w_wid = h5["w_wid"][:]
    d_wid = h5["d_wid"][:]

    line_centre_unit = u"cm^-1"


    #The following is from full_physics/forwardmodel/sunspect
    margin=100

    SOLAR_ANGULAR_RADIUS = (959.44/ 3600) * (pi / 180)

    fovo = 9.2e-3
    frac = fovo / (2 * SOLAR_ANGULAR_RADIUS) #Fraction of the solar diameter viewed, equation from calc_solar


    #solar limb darkening?
    #sld=2/(1+sqrt(1-frac^22))
    sld=1 #UoL-FP uses this one

    #Select the solar needed for this spectral window (with a margin)
    kline1=searchsortedfirst(freq,swin.ww_grid[1]-margin)-1
    kline2=searchsortedfirst(freq,swin.ww_grid[end]+margin)-1

    transmittance = zeros(swin.N_hires)

    for line = kline1:kline2
        if stren[line] < 0
            this_str= 0
        else
            this_str = stren[line]
        end
        this_freq = freq[line]
        this_w_wid = w_wid[line]
        this_d_wid = d_wid[line]
        srot=5e-06*this_freq*frac #broadening due to solar rotation
        d4=(this_d_wid^2+srot^2)^2  # Total Gaussian width
        flinwid=sqrt(2*this_str*(this_d_wid+this_w_wid)/0.0001)

        #if the broadened line lies outside our wavenumber range, we don't need to consider it
        if ((this_freq + flinwid) < swin.ww_grid[1])  continue end
        if ((this_freq - flinwid) > swin.ww_grid[end])  continue end

        y2=(this_w_wid)^2
        ss=sld*this_str
        for iv= 1:swin.N_hires
            xx = swin.ww_grid[iv] - this_freq
            if (abs(xx) > flinwid) continue end
            x2=xx^2
            rr=x2/sqrt(d4+y2*x2*(1+abs(xx/(this_w_wid+0.07))))
            yy=ss*exp(-rr)
            transmittance[iv]-=yy
        end
    end

    transmittance = exp.(transmittance)

    #The following is from full_physics/forwardmodel/calc_solar

    #continuum spectrum is calculated in ph/s/m2/micron so we need to convert wavenumber to microns
    cont_microns=1e4./swin.ww_grid

    bb = [-7.0251527e+22,3.1243395e+23,-4.2464027e+23,1.8903014e+23] #these values are from /data/ghgas/GOSAT/input/template/static_input/in/solar/solar_v2_2019.dat

    continuum = zeros(length(cont_microns))
    for (i,this_wl) in enumerate(cont_microns)
        for j = 1:4
            continuum[i]+=bb[j]*this_wl^(j-1)
        end
    end

    #continuum radiance is in units ph/s/m2/micron so convert to W/m2/cm-1

    #1) convert ph/s into W
    @views continuum[:] .*= ustrip.(Ref(u"W"),
        1.0u"s^-1" .* SPEED_OF_LIGHT ./ (cont_microns .* u"µm") .* PLANCK
    )

    #2) convert  W/m2/µm into W/m2/cm^-1
    @views continuum[:] ./= (1e4 ./ cont_microns) .^ 2

    continuum = continuum[end:-1:1]

    ww_unit = u"cm^-1"
    ww = swin.ww_grid

    irradiance_unit = u"W/m^2/cm^-1" 

    # Close up HDF file, all done
    close(h5)

    # Return solar model object
    return UoLFPSolarModel(
        filename,
        ww,
        transmittance,
        continuum,
        ww_unit,
        irradiance_unit
    )

end


"""
Calculates the solar Doppler shift factor for a given
`EarthLocation` and `DateTime`.

$(TYPEDSIGNATURES)

"""
function calculate_solar_doppler_shift(
    scene::EarthScene
    )

    #Calculate portion of doppler shift due to rotation of earth away from sun
    earth_radius = 6378.137e3
    earth_rot_freq = 2*pi/86164.09054#earth's angular rotation frequency

    geocen_lat = atand(tand(scene.location.latitude)/(1+6.73951496e-3)) #calculate geocentric latitude from geodetic latitude
    gcrad = scene.location.altitude + earth_radius/sqrt(1+6.73951496e-3*sind(geocen_lat)^2) #local earth radius
    earth_rot_velocity = -earth_rot_freq * gcrad * sind(scene.solar_zenith) * cosd(scene.solar_azimuth-90u"°") * cosd(geocen_lat)

    #Calculate portion of doppler shift due to movement of earth center away from sun
    a = [-1.82823e-5, 2.30179e-6, 6.62402e-9, -1.33287e-10, 3.98445e-13, -3.54239e-16]
    days_of_year = Dates.toms(scene.time-DateTime(string(year(scene.time)),dateformat"y")) / 1000 /86400 + 1.5
    j = [1,2*days_of_year,3*days_of_year^2,4*days_of_year^3,5*days_of_year^4,6*days_of_year^5] #we want to calculate the derivate wrt time of the distance function
    earth_sun_velocity=(transpose(a)*j)/86400*1.49597870691e11  #convert from AU/day to m/s

    doppler_shift = (earth_rot_velocity + earth_sun_velocity)/2.99792458e8

    return doppler_shift

end

"""
Calculates the Earth-Sun distances for given
 `DateTime`.

$(TYPEDSIGNATURES)

"""

function calculate_earth_sun_distance(time::DateTime)

    #Calculate Earth-Sun distance (in AU)

    #Parameters from 6th order polynomial fit to data from http://eclipse.gsfc.nasa.gov/TYPE/TYPE.html
    a = [0.98334, -1.82823e-5, 2.30179e-6, 6.62402e-9, -1.33287e-10, 3.98445e-13, -3.54239e-16]

    day_of_year = Dates.toms(time-DateTime(string(year(time)),dateformat"y")) ./ 1000 ./86400
    j = [1,day_of_year,day_of_year^2,day_of_year^3,day_of_year^4,day_of_year^5,day_of_year^6]
    solar_earth_dist=(transpose(a)*j)
    return solar_earth_dist
end

#function calculate_earth_sun_distance(time::Vector{DateTime})

    #Calculate Earth-Sun distance (in AU)

    #Parameters from 6th order polynomial fit to data from http://eclipse.gsfc.nasa.gov/TYPE/TYPE.html
#    a = [0.98334, -1.82823e-5, 2.30179e-6, 6.62402e-9, -1.33287e-10, 3.98445e-13, -3.54239e-16]

#    day_of_year = Dates.toms.(time-DateTime.(string.(year.(time)),dateformat"y")) ./ 1000 ./86400
#    j = [ones(length(day_of_year)),day_of_year,day_of_year.^2,day_of_year.^3,day_of_year.^4,day_of_year.^5,day_of_year.^6]
#    day_of_year_powers=[x[i] for x in values(j), i=1:length(first(j))]
#    solar_earth_dist=(transpose(a)*day_of_year_powers)[1,:]

#    return solar_earth_dist
#end