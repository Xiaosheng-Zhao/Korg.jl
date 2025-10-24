"""
# example_single_spectrum.jl

Simple example demonstrating how to generate a single stellar spectrum
using the Kurucz-Korg pipeline.

This is a minimal working example to get you started quickly.

Author: Generated for JWST M31 spectral analysis
Date: 2025
"""

using Korg
using Plots  # Optional, for visualization

# Include the bridge module
include("KuruczKorgBridge.jl")
using .KuruczKorgBridge

# ============================================================================
# Configuration
# ============================================================================

# Path to your Kurucz model
# Update this path to point to your kurucz1 model weights
KURUCZ_MODEL_PATH = get(ENV, "KURUCZ_MODEL_PATH",
                        "../kurucz1/model/a_one_weights.pt")

# Stellar parameters
Teff = 5777.0      # Effective temperature (K) - Solar value
logg = 4.44        # Surface gravity (log10 cm/s²) - Solar value
Fe_H = 0.0         # Metallicity [Fe/H] - Solar value
alpha_Fe = 0.0     # Alpha enhancement [α/Fe] - Solar value

# Wavelength range (Angstroms)
# Example: JWST NIRSpec range suitable for M31 stellar spectroscopy
wavelengths = [(15000.0, 28000.0, 0.5)]  # 1.5-2.8 μm with 0.5 Å spacing

# Microturbulence velocity (km/s)
vmic = 1.0

# ============================================================================
# Step-by-step workflow
# ============================================================================

println("="^70)
println("Single Spectrum Generation Example")
println("="^70)

# Step 1: Load the Kurucz emulator
println("\n[1/4] Loading Kurucz atmosphere emulator...")
emulator = load_kurucz_emulator(KURUCZ_MODEL_PATH)
println("✓ Emulator loaded successfully")

# Step 2: Generate Kurucz atmospheric structure
println("\n[2/4] Generating Kurucz atmospheric structure...")
println("  Parameters: Teff=$Teff K, logg=$logg, [Fe/H]=$Fe_H, [α/Fe]=$alpha_Fe")
kurucz_atm = KuruczKorgBridge.get_kurucz_atmosphere(
    emulator, Teff, logg, Fe_H, alpha_Fe; n_tau=80
)
println("✓ Atmospheric structure generated ($(length(kurucz_atm.tau)) layers)")

# Step 3: Convert to Korg format
println("\n[3/4] Converting to Korg atmosphere format...")
korg_atm = kurucz_to_korg(kurucz_atm)
println("✓ Conversion complete")

# Step 4: Synthesize spectrum
println("\n[4/4] Synthesizing stellar spectrum...")
println("  Wavelength range: $(wavelengths[1][1])-$(wavelengths[1][2]) Å")
println("  Microturbulence: $vmic km/s")

# Prepare abundances
A_X = Korg.format_A_X(Fe_H, alpha_H=alpha_Fe)

# Synthesize
result = synthesize_spectrum(korg_atm, A_X, wavelengths; vmic=vmic)

println("✓ Synthesis complete!")
println("  Generated $(length(result.wavelengths)) wavelength points")

# ============================================================================
# Results and output
# ============================================================================

println("\n" * "="^70)
println("Results Summary")
println("="^70)
println("Wavelength range: $(minimum(result.wavelengths)) - $(maximum(result.wavelengths)) Å")
println("Number of points: $(length(result.wavelengths))")
println("Flux range: $(minimum(result.flux)) - $(maximum(result.flux)) erg/s/cm²/Å")
println()

# Optional: Save to file
output_file = "spectrum_$(Int(Teff))K_logg$(logg)_FeH$(Fe_H).dat"
open(output_file, "w") do io
    println(io, "# Stellar spectrum generated with Kurucz-Korg pipeline")
    println(io, "# Teff=$Teff K, logg=$logg, [Fe/H]=$Fe_H, [α/Fe]=$alpha_Fe")
    println(io, "# Wavelength(Angstrom)  Flux  Continuum")
    for i in 1:length(result.wavelengths)
        println(io, "$(result.wavelengths[i])  $(result.flux[i])  $(result.continuum[i])")
    end
end
println("Spectrum saved to: $output_file")

# Optional: Quick plot
try
    using Plots
    plt = plot(result.wavelengths, result.flux,
               xlabel="Wavelength (Å)",
               ylabel="Flux",
               title="Synthetic Spectrum: Teff=$Teff K, [Fe/H]=$Fe_H",
               label="Spectrum",
               linewidth=1)
    plot!(plt, result.wavelengths, result.continuum,
          label="Continuum", linewidth=2, linestyle=:dash)
    savefig(plt, "spectrum_plot.png")
    println("Plot saved to: spectrum_plot.png")
catch
    println("Note: Install Plots.jl to generate visualization (optional)")
end

println("\n" * "="^70)
println("Example complete!")
println("="^70)

# ============================================================================
# Alternative: One-line convenience function
# ============================================================================

println("\n\nAlternative: Using the convenience function...")
println("="^70)

# You can also do everything in one function call:
result2 = KuruczKorgBridge.kurucz_synthesize(
    emulator, Teff, logg, Fe_H, alpha_Fe, wavelengths; vmic=vmic
)

println("✓ One-line synthesis complete!")
println("  Same results: $(result.wavelengths ≈ result2.wavelengths)")
