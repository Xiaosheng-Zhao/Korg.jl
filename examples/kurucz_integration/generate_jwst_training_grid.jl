"""
# generate_jwst_training_grid.jl

Generate a grid of synthetic stellar spectra for JWST wavelength ranges.
This creates training data for stellar parameter inference from JWST M31 observations.

Usage:
    julia generate_jwst_training_grid.jl

This script will:
1. Load the Kurucz atmosphere emulator
2. Generate a grid of stellar parameters (Teff, logg, [Fe/H], [α/Fe])
3. For each parameter set, synthesize spectra on JWST wavelength grids
4. Save results to HDF5 or FITS files for machine learning training

Author: Generated for JWST M31 spectral analysis
Date: 2025
"""

using Korg
using HDF5
using ProgressMeter
using Dates

# Include the bridge module
include("KuruczKorgBridge.jl")
using .KuruczKorgBridge

# ============================================================================
# JWST Wavelength Configurations
# ============================================================================

"""
JWST NIRSpec wavelength ranges (in Angstroms)

References:
- NIRSpec Prism: 6000-53000 Å (0.6-5.3 μm)
- NIRSpec G140M: 9700-18400 Å (0.97-1.84 μm)
- NIRSpec G235M: 17000-30700 Å (1.70-3.07 μm)
- NIRSpec G395M: 29000-51800 Å (2.90-5.18 μm)
- NIRSpec G140H: 9700-18400 Å (high resolution)
- NIRSpec G235H: 17000-30700 Å (high resolution)
- NIRSpec G395H: 29000-51800 Å (high resolution)
"""
const JWST_WAVELENGTH_CONFIGS = Dict(
    # Full optical to near-IR range (for reference, outside JWST)
    :optical => [(4000.0, 10000.0, 0.1)],

    # JWST NIRSpec configurations
    :nirspec_prism => [(6000.0, 53000.0, 5.0)],  # Lower resolution for speed
    :nirspec_g140m => [(9700.0, 18400.0, 0.5)],
    :nirspec_g235m => [(17000.0, 30700.0, 0.5)],
    :nirspec_g395m => [(29000.0, 51800.0, 0.5)],
    :nirspec_full => [(9700.0, 51800.0, 1.0)],   # Combined NIRSpec range

    # JWST NIRCam (for photometric cross-calibration)
    :nircam_sw => [(6000.0, 23500.0, 2.0)],  # Short wavelength
    :nircam_lw => [(23500.0, 50000.0, 2.0)], # Long wavelength

    # Custom M31 survey range (adjust based on your specific needs)
    :m31_survey => [(15000.0, 28000.0, 0.5)],  # Common range for stellar spectroscopy
)

# ============================================================================
# Stellar Parameter Grid Definition
# ============================================================================

"""
Define stellar parameter grid for M31 stellar populations.

Grid ranges chosen to cover:
- Main sequence dwarfs: High Teff, high logg
- Red giants: Lower Teff, low logg
- Various metallicities: Thin disk, thick disk, halo populations
- Alpha enhancements: Different star formation histories

Adjust these ranges based on your specific science goals.
"""
function get_stellar_parameter_grid(;
    # Temperature range (K)
    Teff_min=3500.0, Teff_max=7000.0, Teff_step=250.0,

    # Surface gravity (log10(cm/s²))
    logg_min=1.0, logg_max=5.0, logg_step=0.5,

    # Metallicity [Fe/H]
    FeH_min=-2.5, FeH_max=0.5, FeH_step=0.5,

    # Alpha enhancement [α/Fe]
    alphaFe_min=0.0, alphaFe_max=0.4, alphaFe_step=0.2
)

    # Create grid
    Teff_grid = Teff_min:Teff_step:Teff_max
    logg_grid = logg_min:logg_step:logg_max
    FeH_grid = FeH_min:FeH_step:FeH_max
    alphaFe_grid = alphaFe_min:alphaFe_step:alphaFe_max

    # Generate all combinations
    grid = []
    for Teff in Teff_grid
        for logg in logg_grid
            for FeH in FeH_grid
                for alphaFe in alphaFe_grid
                    push!(grid, (Teff=Teff, logg=logg, FeH=FeH, alphaFe=alphaFe))
                end
            end
        end
    end

    @info "Stellar parameter grid created" n_models=length(grid) n_Teff=length(Teff_grid) n_logg=length(logg_grid) n_FeH=length(FeH_grid) n_alphaFe=length(alphaFe_grid)

    return grid
end

"""
Get a smaller test grid for quick validation
"""
function get_test_grid()
    return [
        (Teff=5777.0, logg=4.44, FeH=0.0, alphaFe=0.0),    # Sun
        (Teff=5000.0, logg=4.5, FeH=-0.5, alphaFe=0.2),    # Metal-poor dwarf
        (Teff=4500.0, logg=2.5, FeH=-1.0, alphaFe=0.3),    # Metal-poor giant
        (Teff=6000.0, logg=4.0, FeH=0.0, alphaFe=0.0),     # Hot dwarf
        (Teff=4000.0, logg=1.5, FeH=-0.5, alphaFe=0.2),    # Cool giant
    ]
end

# ============================================================================
# Spectrum Generation Functions
# ============================================================================

"""
    generate_spectrum_for_params(emulator, params, wavelengths; kwargs...)

Generate spectrum for a single set of stellar parameters.

# Arguments
- `emulator`: Loaded Kurucz emulator
- `params::NamedTuple`: (Teff, logg, FeH, alphaFe)
- `wavelengths`: Wavelength configuration
- `kwargs...`: Additional arguments (vmic, linelist, etc.)

# Returns
- `NamedTuple`: (params=..., wavelengths=..., flux=..., continuum=...)
"""
function generate_spectrum_for_params(emulator, params, wavelengths; vmic=1.0, linelist=nothing)
    try
        result = KuruczKorgBridge.kurucz_synthesize(
            emulator,
            params.Teff,
            params.logg,
            params.FeH,
            params.alphaFe,
            wavelengths;
            vmic=vmic,
            linelist=linelist
        )

        return (
            params = params,
            wavelengths = result.wavelengths,
            flux = result.flux,
            continuum = result.continuum,
            success = true,
            error_msg = ""
        )
    catch e
        @warn "Failed to generate spectrum" params error=e
        return (
            params = params,
            wavelengths = nothing,
            flux = nothing,
            continuum = nothing,
            success = false,
            error_msg = string(e)
        )
    end
end

"""
    generate_grid(emulator, param_grid, wavelengths; output_file, kwargs...)

Generate spectra for entire parameter grid and save to HDF5 file.

# Arguments
- `emulator`: Loaded Kurucz emulator
- `param_grid::Vector`: List of parameter sets
- `wavelengths`: Wavelength configuration
- `output_file::String`: Output HDF5 file path
- `vmic::Real`: Microturbulence velocity (km/s, default: 1.0)
- `linelist`: Atomic/molecular line list (optional)
- `save_failed::Bool`: Save failed synthesis attempts (default: false)

# Returns
- Nothing (saves to HDF5 file)
"""
function generate_grid(emulator, param_grid, wavelengths;
                      output_file::String="jwst_training_spectra.h5",
                      vmic::Real=1.0,
                      linelist=nothing,
                      save_failed::Bool=false)

    @info "Starting grid generation" n_models=length(param_grid) output_file

    # Load linelist once (more efficient)
    if linelist === nothing
        @info "Loading VALD solar linelist..."
        linelist = Korg.get_VALD_solar_linelist()
    end

    # Storage for results
    successful_results = []
    failed_results = []

    # Progress bar
    p = Progress(length(param_grid), desc="Generating spectra: ")

    # Generate spectra
    for params in param_grid
        result = generate_spectrum_for_params(emulator, params, wavelengths;
                                             vmic=vmic, linelist=linelist)

        if result.success
            push!(successful_results, result)
        else
            push!(failed_results, result)
        end

        next!(p)
    end

    # Report statistics
    n_success = length(successful_results)
    n_failed = length(failed_results)
    @info "Grid generation complete" n_success n_failed success_rate=n_success/(n_success+n_failed)

    # Save to HDF5
    if n_success > 0
        save_to_hdf5(successful_results, output_file)
        @info "Saved successful results to $output_file"
    end

    if save_failed && n_failed > 0
        failed_file = replace(output_file, ".h5" => "_failed.txt")
        save_failed_params(failed_results, failed_file)
        @info "Saved failed parameters to $failed_file"
    end

    return (successful=successful_results, failed=failed_results)
end

"""
    save_to_hdf5(results, filename)

Save spectral grid to HDF5 file for machine learning training.

HDF5 structure:
- /wavelengths: Common wavelength grid [n_wavelength]
- /spectra: Flux array [n_models, n_wavelength]
- /continua: Continuum array [n_models, n_wavelength]
- /parameters/Teff: [n_models]
- /parameters/logg: [n_models]
- /parameters/FeH: [n_models]
- /parameters/alphaFe: [n_models]
- /metadata: Creation date, code version, etc.
"""
function save_to_hdf5(results, filename::String)
    n_models = length(results)

    # Assume all results have the same wavelength grid
    wavelengths = results[1].wavelengths

    # Preallocate arrays
    n_wl = length(wavelengths)
    spectra = zeros(Float64, n_models, n_wl)
    continua = zeros(Float64, n_models, n_wl)
    Teffs = zeros(Float64, n_models)
    loggs = zeros(Float64, n_models)
    FeHs = zeros(Float64, n_models)
    alphaFes = zeros(Float64, n_models)

    # Fill arrays
    for (i, result) in enumerate(results)
        spectra[i, :] = result.flux
        continua[i, :] = result.continuum
        Teffs[i] = result.params.Teff
        loggs[i] = result.params.logg
        FeHs[i] = result.params.FeH
        alphaFes[i] = result.params.alphaFe
    end

    # Write to HDF5
    h5open(filename, "w") do file
        # Wavelengths
        file["wavelengths"] = wavelengths
        attrs(file["wavelengths"])["units"] = "Angstrom"

        # Spectra
        file["spectra"] = spectra
        attrs(file["spectra"])["units"] = "erg/s/cm^2/cm (specific intensity)"
        attrs(file["spectra"])["description"] = "Synthetic stellar spectra [n_models, n_wavelength]"

        file["continua"] = continua
        attrs(file["continua"])["units"] = "erg/s/cm^2/cm (specific intensity)"

        # Parameters
        g = create_group(file, "parameters")
        g["Teff"] = Teffs
        attrs(g["Teff"])["units"] = "K"
        g["logg"] = loggs
        attrs(g["logg"])["units"] = "log10(cm/s^2)"
        g["FeH"] = FeHs
        attrs(g["FeH"])["units"] = "dex"
        g["alphaFe"] = alphaFes
        attrs(g["alphaFe"])["units"] = "dex"

        # Metadata
        m = create_group(file, "metadata")
        m["creation_date"] = string(now())
        m["n_models"] = n_models
        m["n_wavelength"] = n_wl
        m["generator"] = "Kurucz-Korg pipeline"
        m["atmosphere_model"] = "Kurucz neural emulator"
        m["synthesis_code"] = "Korg.jl"
    end

    @info "Saved $n_models spectra to HDF5" filename size_MB=filesize(filename)/1e6
end

"""
    save_failed_params(failed_results, filename)

Save list of failed parameter combinations to text file.
"""
function save_failed_params(failed_results, filename::String)
    open(filename, "w") do io
        println(io, "# Failed parameter combinations")
        println(io, "# Teff, logg, FeH, alphaFe, Error")
        for result in failed_results
            p = result.params
            println(io, "$(p.Teff), $(p.logg), $(p.FeH), $(p.alphaFe), $(result.error_msg)")
        end
    end
end

# ============================================================================
# Main execution
# ============================================================================

"""
Main function to run the grid generation.

Modify the parameters below to customize your training grid.
"""
function main()
    println("="^70)
    println("JWST M31 Training Data Generation")
    println("Kurucz Atmospheres + Korg Spectral Synthesis")
    println("="^70)

    # ========================================================================
    # CONFIGURATION - Modify these parameters as needed
    # ========================================================================

    # Path to Kurucz model weights
    kurucz_model_path = ENV["KURUCZ_MODEL_PATH"] = get(ENV, "KURUCZ_MODEL_PATH",
                                                        "../kurucz1/model/a_one_weights.pt")

    # Output file
    output_file = "jwst_m31_training_spectra.h5"

    # Choose wavelength configuration
    wavelength_config = :m31_survey  # Options: :nirspec_prism, :nirspec_full, :m31_survey, etc.
    wavelengths = JWST_WAVELENGTH_CONFIGS[wavelength_config]

    # Microturbulence velocity (km/s)
    vmic = 1.0

    # Use test grid or full grid
    use_test_grid = true  # Set to false for full production run

    # ========================================================================
    # EXECUTION
    # ========================================================================

    # Load Kurucz emulator
    @info "Loading Kurucz emulator from $kurucz_model_path"
    emulator = load_kurucz_emulator(kurucz_model_path)

    # Get parameter grid
    if use_test_grid
        @info "Using test grid (5 models)"
        param_grid = get_test_grid()
    else
        @info "Using full parameter grid"
        param_grid = get_stellar_parameter_grid(
            Teff_min=3500.0, Teff_max=7000.0, Teff_step=250.0,
            logg_min=1.0, logg_max=5.0, logg_step=0.5,
            FeH_min=-2.5, FeH_max=0.5, FeH_step=0.5,
            alphaFe_min=0.0, alphaFe_max=0.4, alphaFe_step=0.2
        )
    end

    # Generate grid
    @info "Wavelength configuration" config=wavelength_config ranges=wavelengths
    results = generate_grid(emulator, param_grid, wavelengths;
                          output_file=output_file,
                          vmic=vmic,
                          save_failed=true)

    # Print summary
    println("\n" * "="^70)
    println("Grid generation complete!")
    println("  Successful: $(length(results.successful))")
    println("  Failed: $(length(results.failed))")
    println("  Output: $output_file")
    println("="^70)

    return results
end

# Run main if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
