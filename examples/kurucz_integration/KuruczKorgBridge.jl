"""
# KuruczKorgBridge.jl

Module to bridge Kurucz stellar atmosphere emulator (Python) with Korg.jl (Julia)
for stellar spectra generation.

This module provides functions to:
1. Interface with the Python-based Kurucz emulator
2. Convert Kurucz atmospheric structure to Korg format
3. Generate stellar spectra using Korg synthesis

Author: Generated for JWST M31 spectral analysis
Date: 2025
"""

module KuruczKorgBridge

using PyCall
using Korg

# ----------------------------------------------------------------------
# Universal helper: works whether or not pyconvert/py_to_array exist
# ----------------------------------------------------------------------
function py_to_array(x)
    # Case 1: Already a Julia Array
    if isa(x, AbstractArray)
        return x
    end

    # Case 2: PyTorch tensor (has .detach or .numpy)
    if pyhasattr(x, "detach")
        return x.detach().cpu().numpy()
    elseif pyhasattr(x, "numpy")
        return x.numpy()
    end

    # Case 3: NumPy array (we just want to turn it into a Julia Array)
    try
        np = pyimport("numpy")
        return Array(np.asarray(x))
    catch
        # Worst-case fallback
        return Array(x)
    end
end


export KuruczAtmosphere, kurucz_to_korg, synthesize_spectrum, load_kurucz_emulator

# Physical constants
const k_B = 1.380649e-16  # Boltzmann constant (erg/K)
const m_H = 1.6735575e-24  # Hydrogen atom mass (g)
const μ = 1.3  # Mean molecular weight (approximate for stellar atmospheres)

"""
    load_kurucz_emulator(model_path, norm_params_path=nothing)

Load the Kurucz atmosphere emulator using PyCall.

# Arguments
- `model_path::String`: Path to the Kurucz model weights (e.g., "a_one_weights.pt")
- `norm_params_path::Union{String,Nothing}`: Path to normalization parameters (optional)

# Returns
- Python emulator object

# Example
```julia
emulator = load_kurucz_emulator("../model/a_one_weights.pt")
```
"""
function load_kurucz_emulator(model_path::String, norm_params_path=nothing)
    try
        kuruczone = pyimport("kuruczone")
        emulator_module = kuruczone.emulator

        if norm_params_path === nothing
            model = emulator_module.load_from_checkpoint(model_path)
        else
            model = emulator_module.load_from_checkpoint(model_path, norm_params_path)
        end

        @info "Kurucz emulator loaded successfully from $model_path"
        return model
    catch e
        error("""
        Failed to load Kurucz emulator. Please ensure:
        1. Python environment has kuruczone installed: pip install -e path/to/kurucz1
        2. PyTorch is installed: pip install torch
        3. The model path is correct: $model_path

        Error: $e
        """)
    end
end

"""
    KuruczAtmosphere

Struct to hold Kurucz atmospheric structure data.

# Fields
- `stellar_params::Vector{Float64}`: [Teff, logg, [Fe/H], [α/Fe]]
- `tau::Vector{Float64}`: Optical depth grid
- `T::Vector{Float64}`: Temperature (K)
- `P::Vector{Float64}`: Gas pressure (dyn/cm²)
- `RHOX::Vector{Float64}`: Mass column density (g/cm²)
- `XNE::Vector{Float64}`: Electron number density (cm⁻³)
- `ABROSS::Vector{Float64}`: Rosseland mean absorption coefficient (cm²/g)
- `ACCRAD::Vector{Float64}`: Radiative acceleration (cm/s²)
"""
struct KuruczAtmosphere
    stellar_params::Vector{Float64}
    tau::Vector{Float64}
    T::Vector{Float64}
    P::Vector{Float64}
    RHOX::Vector{Float64}
    XNE::Vector{Float64}
    ABROSS::Vector{Float64}
    ACCRAD::Vector{Float64}
end

"""
    get_kurucz_atmosphere(emulator, Teff, logg, Fe_H, alpha_Fe; n_tau=80)

Generate atmospheric structure using Kurucz emulator.

# Arguments
- `emulator`: Loaded Kurucz emulator object
- `Teff::Real`: Effective temperature (K)
- `logg::Real`: Surface gravity (log10(cm/s²))
- `Fe_H::Real`: Metallicity [Fe/H]
- `alpha_Fe::Real`: Alpha enhancement [α/Fe]
- `n_tau::Int`: Number of optical depth points (default: 80)

# Returns
- `KuruczAtmosphere`: Atmospheric structure

# Example
```julia
emulator = load_kurucz_emulator("model.pt")
kurucz_atm = get_kurucz_atmosphere(emulator, 5777.0, 4.44, 0.0, 0.0)
```
"""
function get_kurucz_atmosphere(emulator, Teff::Real, logg::Real, Fe_H::Real, alpha_Fe::Real;
                                n_tau::Int=80)
    torch = pyimport("torch")

    # Prepare stellar parameters (as float64, though model converts internally)
    stellar_params = torch.tensor([[Float64(Teff), Float64(logg), Float64(Fe_H), Float64(alpha_Fe)]])

    # Create optical depth grid (log-spaced from 10^-7 to 10^2)
    tau_grid = torch.logspace(-7, 2, n_tau).unsqueeze(0)

    # Predict atmospheric structure (returns NumPy arrays)
    atmosphere = emulator.predict(stellar_params, tau_grid)
    tau     = vec(py_to_array(atmosphere["TAU"]))
    T       = vec(py_to_array(atmosphere["T"]))
    P       = vec(py_to_array(atmosphere["P"]))
    RHOX    = vec(py_to_array(atmosphere["RHOX"]))
    XNE     = vec(py_to_array(atmosphere["XNE"]))
    ABROSS  = vec(py_to_array(atmosphere["ABROSS"]))
    ACCRAD  = vec(py_to_array(atmosphere["ACCRAD"]))

    return KuruczAtmosphere(
        [Teff, logg, Fe_H, alpha_Fe],
        tau, T, P, RHOX, XNE, ABROSS, ACCRAD
    )
end

"""
    calculate_heights(P, RHOX, tau, logg)

Calculate geometric heights from pressure and column mass density.

Uses hydrostatic equilibrium: dP = -g * dρ * dz
Rearranging: dz = -dP / (g * ρ) = -dP / g * (dx/dρ)

# Arguments
- `P::Vector{Float64}`: Gas pressure (dyn/cm²)
- `RHOX::Vector{Float64}`: Mass column density (g/cm²)
- `tau::Vector{Float64}`: Optical depth at reference wavelength
- `logg::Real`: Surface gravity (log10(cm/s²))

# Returns
- `z::Vector{Float64}`: Heights (cm), normalized so z=0 at τ_ross ≈ 1
"""
function calculate_heights(P::Vector{Float64}, RHOX::Vector{Float64},
                          tau::Vector{Float64}, logg::Real)
    g = 10.0^logg  # cm/s²
    n = length(P)
    z = zeros(Float64, n)

    # Integrate from top (low pressure) to bottom (high pressure)
    # z[i] = ∫ (1/ρ) dP/g where ρ = dRHOX/dz
    for i in 2:n
        dRHOX = RHOX[i] - RHOX[i-1]
        dP = P[i] - P[i-1]

        if dRHOX != 0
            # ρ ≈ dRHOX/dz, so dz = dRHOX/ρ
            # From hydrostatic: dP = g * ρ * dz, so dz = dP/(g*ρ)
            # Average density between layers
            rho_avg = abs(dRHOX) / (abs(dP) / g)
            dz = dP / (g * rho_avg)
            z[i] = z[i-1] + dz
        else
            z[i] = z[i-1]
        end
    end

    # Normalize so z=0 at optical depth = 1 (photosphere)
    # Find index closest to τ = 1
    tau_1_idx = argmin(abs.(tau .- 1.0))
    return z .- z[tau_1_idx]
end

"""
    calculate_number_density(P, T)

Calculate total number density from pressure and temperature using ideal gas law.

P = n * k_B * T
n = P / (k_B * T)

# Arguments
- `P::Vector{Float64}`: Gas pressure (dyn/cm²)
- `T::Vector{Float64}`: Temperature (K)

# Returns
- `n::Vector{Float64}`: Total number density (cm⁻³)
"""
function calculate_number_density(P::Vector{Float64}, T::Vector{Float64})
    return P ./ (k_B .* T)
end

"""
    kurucz_to_korg(kurucz_atm::KuruczAtmosphere; reference_wavelength=5e-5)

Convert Kurucz atmospheric structure to Korg PlanarAtmosphere format.

# Arguments
- `kurucz_atm::KuruczAtmosphere`: Kurucz atmospheric structure
- `reference_wavelength::Real`: Reference wavelength in cm (default: 5e-5 cm = 5000 Å)

# Returns
- `Korg.PlanarAtmosphere`: Atmosphere object for Korg synthesis

# Example
```julia
kurucz_atm = get_kurucz_atmosphere(emulator, 5777.0, 4.44, 0.0, 0.0)
korg_atm = kurucz_to_korg(kurucz_atm)
```
"""
function kurucz_to_korg(kurucz_atm::KuruczAtmosphere; reference_wavelength::Real=5e-5)
    Teff, logg, Fe_H, alpha_Fe = kurucz_atm.stellar_params

    # Calculate derived quantities
    z = calculate_heights(kurucz_atm.P, kurucz_atm.RHOX, kurucz_atm.tau, logg)
    number_density = calculate_number_density(kurucz_atm.P, kurucz_atm.T)

    n_layers = length(kurucz_atm.tau)

    # Let Julia infer the correct element type for PlanarAtmosphereLayer
    layers = [Korg.PlanarAtmosphereLayer(
                  kurucz_atm.tau[i],
                  z[i],
                  kurucz_atm.T[i],
                  kurucz_atm.XNE[i],
                  number_density[i]
              ) for i in 1:n_layers]
    
    # Ensure `reference_wavelength` is Float64 (Å to cm conversion if needed)
    korg_atm = Korg.PlanarAtmosphere(collect(layers), Float64(reference_wavelength))

    @info "Converted Kurucz atmosphere to Korg format" Teff logg Fe_H alpha_Fe n_layers

    return korg_atm
end

"""
    synthesize_spectrum(korg_atm, A_X, wavelengths;
                       vmic=1.0, linelist=nothing,
                       save_continuum=true)

Synthesize stellar spectrum using Korg.

# Arguments
- `korg_atm::Korg.PlanarAtmosphere`: Korg atmosphere object
- `A_X`: Abundance vector (use Korg.format_A_X)
- `wavelengths`: Wavelength range(s) - see Korg documentation
- `vmic::Real`: Microturbulence velocity (km/s, default: 1.0)
- `linelist`: Atomic/molecular line list (default: Korg VALD solar linelist)
- `save_continuum::Bool`: Whether to save continuum (default: true)

# Returns
- `NamedTuple`: (wavelengths=..., flux=..., continuum=..., synthesis_result=...)

# Example
```julia
#A_X = Korg.format_A_X(0.0, alpha_elements=0.0)  # Solar composition
A_X = Korg.format_A_X(0.0, alpha_elements=0.0)
wavelengths = [(9000, 18000, 0.1)]  # JWST NIRSpec range in Å
result = synthesize_spectrum(korg_atm, A_X, wavelengths; vmic=1.2)
```
"""
function synthesize_spectrum(korg_atm::Korg.PlanarAtmosphere, A_X, wavelengths;
                            vmic::Real=1.0, linelist=nothing,
                            save_continuum::Bool=true)
    # Load default linelist if not provided
    if linelist === nothing
        #@info "Loading VALD solar linelist..."
        #linelist = Korg.get_VALD_solar_linelist()
	@info "Loading NIR linelist for JWST wavelength range..."
    	linelist = Korg.read_linelist(
        	get(ENV, "KORG_NIR_LINELIST", "/scratch/zxs/scripts/Korg.jl/examples/kurucz_integration/jwst_nir_linelist_900_1800.dat"),
        	format="kurucz"
    	)
    end

    # Perform synthesis
    @info "Synthesizing spectrum..." wavelengths vmic
    result = Korg.synthesize(korg_atm, linelist, A_X, wavelengths; vmic=vmic)

    # Return results
    return (
        wavelengths = result.wavelengths,
        flux = result.flux,
        continuum = save_continuum ? result.cntm : nothing,
        synthesis_result = result
    )
end

"""
    kurucz_synthesize(emulator, Teff, logg, Fe_H, alpha_Fe, wavelengths;
                      vmic=1.0, linelist=nothing, n_tau=80)

End-to-end function: Kurucz atmosphere generation → Korg synthesis.

# Arguments
- `emulator`: Loaded Kurucz emulator
- `Teff::Real`: Effective temperature (K)
- `logg::Real`: Surface gravity (log10(cm/s²))
- `Fe_H::Real`: Metallicity [Fe/H]
- `alpha_Fe::Real`: Alpha enhancement [α/Fe]
- `wavelengths`: Wavelength range(s) for synthesis
- `vmic::Real`: Microturbulence velocity (km/s, default: 1.0)
- `linelist`: Atomic/molecular line list (optional)
- `n_tau::Int`: Number of optical depth points (default: 80)

# Returns
- `NamedTuple`: Synthesis results with wavelengths and flux

# Example
```julia
emulator = load_kurucz_emulator("model.pt")
result = kurucz_synthesize(emulator, 5777, 4.44, 0.0, 0.0, [(9000, 18000, 0.1)])
```
"""
function kurucz_synthesize(emulator, Teff::Real, logg::Real, Fe_H::Real, alpha_Fe::Real,
                          wavelengths; vmic::Real=1.0, linelist=nothing, n_tau::Int=80)
    # Step 1: Get Kurucz atmosphere
    kurucz_atm = get_kurucz_atmosphere(emulator, Teff, logg, Fe_H, alpha_Fe; n_tau=n_tau)

    # Step 2: Convert to Korg format
    korg_atm = kurucz_to_korg(kurucz_atm)

    # Step 3: Prepare abundances
    A_X = Korg.format_A_X(Fe_H, alpha_elements=alpha_Fe)

    # Step 4: Synthesize spectrum
    result = synthesize_spectrum(korg_atm, A_X, wavelengths; vmic=vmic, linelist=linelist)

    return result
end

end # module
