"""
# visualize_results.jl

Utility script to visualize and inspect generated spectral grids.

Usage:
    julia visualize_results.jl [path_to_hdf5_file]

Author: Generated for JWST M31 spectral analysis
Date: 2025
"""

using HDF5
using Plots
using Statistics
using Printf

"""
    load_spectral_grid(filename)

Load spectral grid from HDF5 file.
"""
function load_spectral_grid(filename::String)
    if !isfile(filename)
        error("File not found: $filename")
    end

    data = h5open(filename, "r") do f
        Dict(
            :wavelengths => read(f["wavelengths"]),
            :spectra => read(f["spectra"]),
            :continua => read(f["continua"]),
            :Teff => read(f["parameters/Teff"]),
            :logg => read(f["parameters/logg"]),
            :FeH => read(f["parameters/FeH"]),
            :alphaFe => read(f["parameters/alphaFe"]),
            :n_models => read(f["metadata/n_models"]),
            :creation_date => read(f["metadata/creation_date"])
        )
    end

    @info "Loaded spectral grid" filename n_models=data[:n_models] creation_date=data[:creation_date]

    return data
end

"""
    print_summary(data)

Print summary statistics of the spectral grid.
"""
function print_summary(data::Dict)
    println("\n" * "="^70)
    println("SPECTRAL GRID SUMMARY")
    println("="^70)

    println("\nDimensions:")
    println("  Models: $(data[:n_models])")
    println("  Wavelength points: $(length(data[:wavelengths]))")
    println("  Wavelength range: $(@sprintf("%.1f", minimum(data[:wavelengths]))) - $(@sprintf("%.1f", maximum(data[:wavelengths]))) Å")
    println("  Wavelength range: $(@sprintf("%.2f", minimum(data[:wavelengths])/1e4)) - $(@sprintf("%.2f", maximum(data[:wavelengths])/1e4)) μm")

    println("\nParameter Ranges:")
    println("  Teff: $(@sprintf("%.0f", minimum(data[:Teff]))) - $(@sprintf("%.0f", maximum(data[:Teff]))) K")
    println("  logg: $(@sprintf("%.2f", minimum(data[:logg]))) - $(@sprintf("%.2f", maximum(data[:logg])))")
    println("  [Fe/H]: $(@sprintf("%.2f", minimum(data[:FeH]))) - $(@sprintf("%.2f", maximum(data[:FeH])))")
    println("  [α/Fe]: $(@sprintf("%.2f", minimum(data[:alphaFe]))) - $(@sprintf("%.2f", maximum(data[:alphaFe])))")

    println("\nUnique Parameter Values:")
    println("  Teff: $(sort(unique(data[:Teff])))")
    println("  logg: $(sort(unique(data[:logg])))")
    println("  [Fe/H]: $(sort(unique(data[:FeH])))")
    println("  [α/Fe]: $(sort(unique(data[:alphaFe])))")

    println("\nFlux Statistics:")
    println("  Mean flux: $(@sprintf("%.2e", mean(data[:spectra])))")
    println("  Flux range: $(@sprintf("%.2e", minimum(data[:spectra]))) - $(@sprintf("%.2e", maximum(data[:spectra])))")

    println("\nFile Information:")
    println("  Creation date: $(data[:creation_date])")

    println("="^70 * "\n")
end

"""
    plot_sample_spectra(data; n_samples=5)

Plot a random sample of spectra from the grid.
"""
function plot_sample_spectra(data::Dict; n_samples::Int=5)
    n_samples = min(n_samples, data[:n_models])
    indices = rand(1:data[:n_models], n_samples)

    plt = plot(
        xlabel="Wavelength (Å)",
        ylabel="Flux (normalized)",
        title="Sample Spectra from Grid",
        legend=:outertopright,
        size=(1200, 600),
        margin=5Plots.mm
    )

    for idx in indices
        # Normalize spectrum
        spec = data[:spectra][idx, :]
        cont = data[:continua][idx, :]
        normalized = spec ./ cont

        label = @sprintf("Teff=%dK, logg=%.1f, [Fe/H]=%.1f",
                        data[:Teff][idx], data[:logg][idx], data[:FeH][idx])

        plot!(plt, data[:wavelengths], normalized, label=label, alpha=0.7, linewidth=1)
    end

    hline!(plt, [1.0], color=:black, linestyle=:dash, label="Continuum", linewidth=1)

    return plt
end

"""
    plot_teff_sequence(data; logg_target=4.5, FeH_target=0.0)

Plot spectra for a sequence of temperatures at fixed logg and [Fe/H].
"""
function plot_teff_sequence(data::Dict; logg_target::Real=4.5, FeH_target::Real=0.0,
                            alphaFe_target::Real=0.0)
    # Find matching models
    tol = 0.1
    mask = (abs.(data[:logg] .- logg_target) .< tol) .&
           (abs.(data[:FeH] .- FeH_target) .< tol) .&
           (abs.(data[:alphaFe] .- alphaFe_target) .< tol)

    indices = findall(mask)

    if isempty(indices)
        @warn "No models found for logg=$logg_target, [Fe/H]=$FeH_target"
        return nothing
    end

    # Sort by Teff
    Teffs = data[:Teff][indices]
    sort_idx = sortperm(Teffs)
    indices = indices[sort_idx]

    plt = plot(
        xlabel="Wavelength (Å)",
        ylabel="Flux (normalized)",
        title=@sprintf("Temperature Sequence (logg=%.1f, [Fe/H]=%.1f)", logg_target, FeH_target),
        legend=:outertopright,
        size=(1200, 600),
        margin=5Plots.mm
    )

    for idx in indices
        spec = data[:spectra][idx, :]
        cont = data[:continua][idx, :]
        normalized = spec ./ cont

        label = @sprintf("Teff=%dK", data[:Teff][idx])
        plot!(plt, data[:wavelengths], normalized, label=label, alpha=0.7, linewidth=1)
    end

    hline!(plt, [1.0], color=:black, linestyle=:dash, label="Continuum", linewidth=1)

    return plt
end

"""
    plot_metallicity_sequence(data; Teff_target=5000.0, logg_target=4.0)

Plot spectra for a sequence of metallicities at fixed Teff and logg.
"""
function plot_metallicity_sequence(data::Dict; Teff_target::Real=5000.0,
                                   logg_target::Real=4.0, alphaFe_target::Real=0.0)
    # Find matching models
    tol = 100.0  # Teff tolerance
    logg_tol = 0.1
    mask = (abs.(data[:Teff] .- Teff_target) .< tol) .&
           (abs.(data[:logg] .- logg_target) .< logg_tol) .&
           (abs.(data[:alphaFe] .- alphaFe_target) .< logg_tol)

    indices = findall(mask)

    if isempty(indices)
        @warn "No models found for Teff=$Teff_target, logg=$logg_target"
        return nothing
    end

    # Sort by [Fe/H]
    FeHs = data[:FeH][indices]
    sort_idx = sortperm(FeHs)
    indices = indices[sort_idx]

    plt = plot(
        xlabel="Wavelength (Å)",
        ylabel="Flux (normalized)",
        title=@sprintf("Metallicity Sequence (Teff=%dK, logg=%.1f)", Teff_target, logg_target),
        legend=:outertopright,
        size=(1200, 600),
        margin=5Plots.mm
    )

    for idx in indices
        spec = data[:spectra][idx, :]
        cont = data[:continua][idx, :]
        normalized = spec ./ cont

        label = @sprintf("[Fe/H]=%.1f", data[:FeH][idx])
        plot!(plt, data[:wavelengths], normalized, label=label, alpha=0.7, linewidth=1)
    end

    hline!(plt, [1.0], color=:black, linestyle=:dash, label="Continuum", linewidth=1)

    return plt
end

"""
    plot_parameter_distribution(data)

Plot distribution of parameters in the grid.
"""
function plot_parameter_distribution(data::Dict)
    p1 = histogram(data[:Teff], bins=20, xlabel="Teff (K)", ylabel="Count",
                   title="Temperature Distribution", legend=false)

    p2 = histogram(data[:logg], bins=20, xlabel="log(g)", ylabel="Count",
                   title="Gravity Distribution", legend=false)

    p3 = histogram(data[:FeH], bins=20, xlabel="[Fe/H]", ylabel="Count",
                   title="Metallicity Distribution", legend=false)

    p4 = histogram(data[:alphaFe], bins=20, xlabel="[α/Fe]", ylabel="Count",
                   title="Alpha Enhancement Distribution", legend=false)

    plt = plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 800), margin=5Plots.mm)

    return plt
end

"""
    export_subset(data, indices, output_file)

Export a subset of the spectral grid to a new HDF5 file.
"""
function export_subset(data::Dict, indices::Vector{Int}, output_file::String)
    h5open(output_file, "w") do f
        f["wavelengths"] = data[:wavelengths]
        f["spectra"] = data[:spectra][indices, :]
        f["continua"] = data[:continua][indices, :]

        g = create_group(f, "parameters")
        g["Teff"] = data[:Teff][indices]
        g["logg"] = data[:logg][indices]
        g["FeH"] = data[:FeH][indices]
        g["alphaFe"] = data[:alphaFe][indices]

        m = create_group(f, "metadata")
        m["n_models"] = length(indices)
        m["creation_date"] = string(Dates.now())
        m["parent_file"] = get(data, :filename, "unknown")
    end

    @info "Exported subset" n_models=length(indices) output_file
end

# ============================================================================
# Main execution
# ============================================================================

function main()
    # Get filename from command line or use default
    filename = length(ARGS) >= 1 ? ARGS[1] : "jwst_m31_training_spectra.h5"

    if !isfile(filename)
        println("Usage: julia visualize_results.jl [path_to_hdf5_file]")
        println("\nDefault file not found: $filename")
        println("Please specify a valid HDF5 file.")
        return
    end

    # Load data
    println("Loading spectral grid from: $filename")
    data = load_spectral_grid(filename)

    # Print summary
    print_summary(data)

    # Generate plots
    println("Generating visualizations...\n")

    # 1. Sample spectra
    println("[1/4] Plotting sample spectra...")
    p1 = plot_sample_spectra(data; n_samples=5)
    savefig(p1, "sample_spectra.png")
    println("  Saved: sample_spectra.png")

    # 2. Temperature sequence
    println("[2/4] Plotting temperature sequence...")
    p2 = plot_teff_sequence(data; logg_target=4.5, FeH_target=0.0)
    if p2 !== nothing
        savefig(p2, "temperature_sequence.png")
        println("  Saved: temperature_sequence.png")
    end

    # 3. Metallicity sequence
    println("[3/4] Plotting metallicity sequence...")
    p3 = plot_metallicity_sequence(data; Teff_target=5000.0, logg_target=4.0)
    if p3 !== nothing
        savefig(p3, "metallicity_sequence.png")
        println("  Saved: metallicity_sequence.png")
    end

    # 4. Parameter distributions
    println("[4/4] Plotting parameter distributions...")
    p4 = plot_parameter_distribution(data)
    savefig(p4, "parameter_distributions.png")
    println("  Saved: parameter_distributions.png")

    println("\n" * "="^70)
    println("Visualization complete!")
    println("="^70)
end

# Run main if executed as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
