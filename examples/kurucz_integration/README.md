# Kurucz-Korg Integration for JWST M31 Spectral Analysis

This directory contains tools to combine **Kurucz stellar atmosphere models** (via neural network emulator) with **Korg.jl spectral synthesis** to generate training data for JWST M31 stellar parameter inference.

## Overview

### Pipeline Architecture

```
Stellar Parameters (Teff, logg, [Fe/H], [α/Fe])
    ↓
Kurucz Neural Emulator (Python) → Atmospheric Structure (T, P, ρ, etc.)
    ↓
Conversion Layer → Korg Atmosphere Format
    ↓
Korg.jl Synthesis → Stellar Spectrum
    ↓
JWST Wavelength Grids → Training Dataset
```

## Table of Contents

1. [Installation & Setup](#installation--setup)
2. [Quick Start](#quick-start)
3. [Detailed Workflow](#detailed-workflow)
4. [JWST Wavelength Configurations](#jwst-wavelength-configurations)
5. [Generating Training Grids](#generating-training-grids)
6. [Limitations & Cautions](#limitations--cautions)
7. [Troubleshooting](#troubleshooting)
8. [References](#references)

---

## Installation & Setup

### Prerequisites

1. **Julia** (≥ 1.6) with Korg.jl installed
2. **Python** (≥ 3.8) with PyTorch and the Kurucz emulator
3. **PyCall.jl** configured to use your Python environment

### Step 1: Install Kurucz Emulator (Python)

```bash
# Clone the Kurucz emulator repository
git clone https://github.com/jiadonglee/kurucz1.git
cd kurucz1

# Install in editable mode
pip install -e .

# Verify installation
python -c "import kuruczone; print('Kurucz emulator installed successfully')"
```

### Step 2: Download Kurucz Model Weights

Download the pre-trained model weights (e.g., `a_one_weights.pt`) from the Kurucz repository or model hosting service and place them in an accessible location.

```bash
# Example location
mkdir -p ~/models/kurucz
# Download or copy a_one_weights.pt to ~/models/kurucz/
```

### Step 3: Configure PyCall.jl

```julia
# In Julia REPL
using Pkg
Pkg.add("PyCall")

# Point PyCall to your Python environment (if needed)
ENV["PYTHON"] = "/path/to/your/python"  # e.g., from `which python`
Pkg.build("PyCall")

# Test PyCall
using PyCall
py"print('PyCall configured successfully')"
```

### Step 4: Install Additional Julia Dependencies

```julia
using Pkg
Pkg.add(["HDF5", "ProgressMeter", "Plots"])  # Plots is optional
```

### Step 5: Set Environment Variable

```bash
# Add to your ~/.bashrc or ~/.zshrc
export KURUCZ_MODEL_PATH="/path/to/a_one_weights.pt"
```

Or set it in Julia:
```julia
ENV["KURUCZ_MODEL_PATH"] = "/path/to/a_one_weights.pt"
```

---

## Quick Start

### Generate a Single Spectrum

```julia
# Navigate to the examples directory
cd("examples/kurucz_integration")

# Run the example
include("example_single_spectrum.jl")
```

This will:
1. Load the Kurucz emulator
2. Generate a solar-like atmosphere (Teff=5777K, logg=4.44, [Fe/H]=0.0)
3. Convert to Korg format
4. Synthesize spectrum in JWST wavelength range (1.5-2.8 μm)
5. Save results to a text file

**Expected output:** A spectrum file and (optionally) a plot.

### Generate a Training Grid

```julia
# Edit configuration in generate_jwst_training_grid.jl first
# Then run:
include("generate_jwst_training_grid.jl")
```

This will generate a grid of spectra and save to HDF5 format for machine learning.

---

## Detailed Workflow

### Step-by-Step Process

#### 1. Load the Kurucz Emulator

```julia
using .KuruczKorgBridge

emulator = load_kurucz_emulator("/path/to/a_one_weights.pt")
```

#### 2. Generate Kurucz Atmospheric Structure

```julia
# Stellar parameters
Teff = 5000.0      # Effective temperature (K)
logg = 4.5         # Surface gravity (log10 cm/s²)
Fe_H = -0.5        # Metallicity [Fe/H]
alpha_Fe = 0.2     # Alpha enhancement [α/Fe]

# Generate atmosphere (80 optical depth layers)
kurucz_atm = KuruczKorgBridge.get_kurucz_atmosphere(
    emulator, Teff, logg, Fe_H, alpha_Fe; n_tau=80
)
```

**Kurucz Output:**
- `tau`: Optical depth grid (10^-7 to 10^2)
- `T`: Temperature (K)
- `P`: Gas pressure (dyn/cm²)
- `RHOX`: Mass column density (g/cm²)
- `XNE`: Electron number density (cm⁻³)
- `ABROSS`: Rosseland mean absorption coefficient
- `ACCRAD`: Radiative acceleration

#### 3. Convert to Korg Format

```julia
korg_atm = kurucz_to_korg(kurucz_atm)
```

**Conversion process:**
- Calculates geometric heights from hydrostatic equilibrium
- Computes total number density from ideal gas law
- Creates `PlanarAtmosphereLayer` objects for each depth point
- Assembles into Korg `PlanarAtmosphere`

#### 4. Synthesize Spectrum

```julia
using Korg

# Set abundances
A_X = Korg.format_A_X(Fe_H, alpha_H=alpha_Fe)

# Define wavelength range (JWST NIRSpec example)
wavelengths = [(15000.0, 28000.0, 0.5)]  # Å, with 0.5 Å spacing

# Synthesize
result = synthesize_spectrum(korg_atm, A_X, wavelengths; vmic=1.0)
```

**Output:**
- `result.wavelengths`: Wavelength grid (Å)
- `result.flux`: Emergent flux (erg/s/cm²/Å)
- `result.continuum`: Continuum level
- `result.synthesis_result`: Full Korg synthesis object

#### 5. (Alternative) One-Line Convenience Function

```julia
result = KuruczKorgBridge.kurucz_synthesize(
    emulator, Teff, logg, Fe_H, alpha_Fe, wavelengths; vmic=1.0
)
```

---

## JWST Wavelength Configurations

Pre-defined wavelength ranges for JWST instruments are available in `generate_jwst_training_grid.jl`:

### NIRSpec Configurations

| Configuration | Wavelength Range (Å) | Resolution | Use Case |
|--------------|---------------------|------------|----------|
| `:nirspec_prism` | 6,000 - 53,000 | Low (R~100) | Full spectral coverage |
| `:nirspec_g140m` | 9,700 - 18,400 | Medium (R~1000) | H-band spectroscopy |
| `:nirspec_g235m` | 17,000 - 30,700 | Medium (R~1000) | K-band spectroscopy |
| `:nirspec_g395m` | 29,000 - 51,800 | Medium (R~1000) | L/M-band spectroscopy |
| `:nirspec_full` | 9,700 - 51,800 | Medium | Combined range |

### Custom M31 Survey Configuration

```julia
:m31_survey => [(15000.0, 28000.0, 0.5)]  # 1.5-2.8 μm, 0.5 Å spacing
```

Adjust this based on your specific observing program!

### Usage

```julia
wavelengths = JWST_WAVELENGTH_CONFIGS[:nirspec_g235m]
```

---

## Generating Training Grids

### Parameter Grid Definition

The grid can be customized in `generate_jwst_training_grid.jl`:

```julia
param_grid = get_stellar_parameter_grid(
    Teff_min=3500.0, Teff_max=7000.0, Teff_step=250.0,
    logg_min=1.0, logg_max=5.0, logg_step=0.5,
    FeH_min=-2.5, FeH_max=0.5, FeH_step=0.5,
    alphaFe_min=0.0, alphaFe_max=0.4, alphaFe_step=0.2
)
```

**Default grid size:** ~2500 models
- Teff: 15 points (3500-7000 K)
- logg: 9 points (1.0-5.0)
- [Fe/H]: 7 points (-2.5 to 0.5)
- [α/Fe]: 3 points (0.0-0.4)

### Stellar Populations Covered

- **Main sequence dwarfs:** High Teff, high logg
- **Red giants:** Lower Teff, low logg (< 3.5)
- **Metallicity range:** Halo (-2.5), thick disk (-1.0), thin disk (0.0)
- **Alpha enhancement:** Different star formation histories

### Running the Grid Generation

```julia
# For testing (5 models)
include("generate_jwst_training_grid.jl")
# Edit: set use_test_grid = true

# For full production run
# Edit: set use_test_grid = false
include("generate_jwst_training_grid.jl")
```

### Output Format (HDF5)

```
jwst_m31_training_spectra.h5
├── wavelengths [n_wavelength]          # Common wavelength grid (Å)
├── spectra [n_models, n_wavelength]    # Flux array (erg/s/cm²/Å)
├── continua [n_models, n_wavelength]   # Continuum array
├── parameters/
│   ├── Teff [n_models]                 # Effective temperature (K)
│   ├── logg [n_models]                 # Surface gravity
│   ├── FeH [n_models]                  # Metallicity [Fe/H]
│   └── alphaFe [n_models]              # Alpha enhancement [α/Fe]
└── metadata/
    ├── creation_date
    ├── n_models
    └── generator = "Kurucz-Korg pipeline"
```

### Loading Training Data (Python)

```python
import h5py
import numpy as np

with h5py.File('jwst_m31_training_spectra.h5', 'r') as f:
    wavelengths = f['wavelengths'][:]
    spectra = f['spectra'][:]
    Teff = f['parameters/Teff'][:]
    logg = f['parameters/logg'][:]
    FeH = f['parameters/FeH'][:]
    alphaFe = f['parameters/alphaFe'][:]

# Use for training neural networks, random forests, etc.
```

---

## Limitations & Cautions

### ⚠️ Critical Limitations

#### 1. **Kurucz Emulator Training Range**

The neural network emulator is only accurate within its training domain:

- **Teff:** Typically 3500-8000 K (check your specific model)
- **logg:** Typically 0.0-5.0
- **[Fe/H]:** Typically -4.0 to +1.0
- **[α/Fe]:** Typically 0.0 to +0.4

**⚠️ Extrapolation beyond these ranges will give unreliable results!**

**Action:** Verify the training range of your specific Kurucz model before defining your parameter grid.

#### 2. **Atmospheric Structure Conversion**

The conversion from Kurucz to Korg format involves several approximations:

**Height calculation:**
- Assumes simple hydrostatic equilibrium integration
- No explicit treatment of turbulent pressure
- Reference point (z=0) is approximate

**Number density:**
- Uses ideal gas law: n = P/(k_B * T)
- Assumes mean molecular weight μ ≈ 1.3
- May not be accurate for very cool stars or unusual compositions

**⚠️ Impact:** Small errors (~few %) in atmospheric structure can propagate to spectrum

**Recommendation:** For critical science cases, validate against MARCS or PHOENIX atmospheres for a subset of your grid.

#### 3. **Optical Depth Grid**

- **Default:** 80 points log-spaced from τ = 10^-7 to 10^2
- **Issue:** May be insufficient for very high spectral resolution
- **Impact:** Can miss fine structure in steep gradients

**Action:** Increase `n_tau` for high-resolution work (e.g., n_tau=100-120 for R>10,000)

#### 4. **Microturbulence**

- **Current:** Single microturbulence value (vmic) for entire atmosphere
- **Reality:** Microturbulence varies with depth
- **Impact:** Line profiles may be slightly incorrect, especially for strong lines

**Typical values:**
- Dwarfs: 1.0-1.5 km/s
- Giants: 1.5-2.5 km/s

#### 5. **Line List Completeness**

The default VALD solar linelist may be incomplete for:
- Cool stars (< 4000 K): Missing molecular lines
- Metal-poor stars ([Fe/H] < -2): Scaled linelist may be inappropriate
- Near-IR: Some atomic/molecular data uncertain

**⚠️ Missing lines → Continuum overestimation**

**Actions:**
- For cool stars: Consider adding molecular linelists (TiO, H₂O, etc.)
- For metal-poor stars: Use dedicated linelists for low-metallicity regimes
- Validate synthetic spectra against observed standards

#### 6. **Computational Performance**

**Grid generation is computationally intensive:**

Example timing (approximate):
- Single spectrum (1.5-2.8 μm, 0.5 Å): ~10-30 seconds
- Grid of 2500 models: ~7-20 hours (single core)

**Recommendations:**
- Use test grid first (`use_test_grid = true`)
- For large grids, parallelize (see "Parallelization" below)
- Consider coarser wavelength sampling for initial tests

#### 7. **Wavelength Coverage**

**Korg's line data quality varies by wavelength:**
- **Best:** Optical (4000-10,000 Å) - well-studied
- **Good:** Near-IR (1.0-2.5 μm) - reasonably complete
- **Uncertain:** Mid-IR (> 3 μm) - many molecular bands poorly known

**For JWST mid-IR (> 3 μm):** Expect larger uncertainties in synthetic spectra.

#### 8. **Stellar Rotation & Macroturbulence**

**Not currently implemented in this pipeline:**
- Rotational broadening (v sin i)
- Macroturbulence
- Instrumental broadening

**Action:** Apply these as post-processing steps if needed for your science case.

#### 9. **Physical Effects Not Included**

The following physical effects are **not modeled**:

- **3D hydrodynamics:** Atmospheres are 1D plane-parallel
- **NLTE:** All calculations assume local thermodynamic equilibrium (LTE)
- **Sphericity (for giants):** Uses plane-parallel geometry (can be severe for log g < 1.5)
- **Magnetic fields**
- **Stellar activity:** Spots, plages, etc.
- **Binary companions**

**Impact:** Systematic errors for:
- Cool giants (3D/NLTE important)
- Strong resonance lines (NLTE important)
- Very metal-poor stars (NLTE important)

#### 10. **M31 vs Milky Way**

**Important considerations for M31 observations:**

- **Abundance patterns:** M31 may have different [X/Fe] ratios than MW
  - **Action:** Adjust abundance patterns using `Korg.format_A_X` if known

- **Interstellar extinction:** Both Galactic and M31 internal
  - **Not included:** Apply extinction corrections separately

- **Radial velocity:** Doppler shifts from M31 rotation
  - **Not included:** Apply wavelength shifts as needed

- **Stellar multiplicity:** Higher in crowded regions
  - **Not included:** Consider composite spectra in data preparation

### 🔍 Validation Recommendations

1. **Compare to Standards:**
   - Generate spectra for well-known stars (e.g., Sun, Arcturus)
   - Compare to observed high-quality spectra
   - Check line depths, continuum shape, features

2. **Cross-Check with Other Codes:**
   - Compare subset of grid to MARCS/PHOENIX + Turbospectrum/MOOG
   - Verify agreement within acceptable tolerances

3. **Parameter Recovery Tests:**
   - Synthesize spectra with known parameters
   - Apply your inference method
   - Check for systematic biases

4. **Uncertainty Quantification:**
   - Propagate uncertainties in atmosphere structure
   - Test sensitivity to linelist changes
   - Assess impact of missing physics (NLTE, 3D, etc.)

---

## Troubleshooting

### Issue: "Failed to load Kurucz emulator"

**Solutions:**
1. Verify Python environment: `python -c "import kuruczone"`
2. Check PyCall configuration: `using PyCall; pyimport("sys")`
3. Rebuild PyCall: `using Pkg; Pkg.build("PyCall")`
4. Check model path: `isfile(KURUCZ_MODEL_PATH)`

### Issue: "Synthesis fails for certain parameter combinations"

**Possible causes:**
1. Parameters outside emulator training range
2. Atmospheric structure is unphysical (check T, P profiles)
3. Insufficient optical depth sampling

**Solutions:**
- Check `failed_params.txt` for patterns
- Narrow parameter range
- Increase `n_tau`
- Validate atmosphere structure before synthesis

### Issue: "Spectra look wrong (too much/too little absorption)"

**Possible causes:**
1. Incorrect abundance scaling
2. Missing molecular lines (cool stars)
3. Wrong microturbulence
4. Height scale calculation error

**Solutions:**
- Compare to reference spectrum (e.g., solar)
- Check `A_X` abundance vector
- Validate atmosphere conversion
- Try different vmic values

### Issue: "Very slow generation"

**Solutions:**
1. Reduce wavelength sampling (larger step size)
2. Reduce spectral range
3. Use test grid first
4. Implement parallelization (see below)

### Parallelization (Advanced)

For large grids, parallelize using Julia's `Distributed`:

```julia
using Distributed
addprocs(8)  # Use 8 cores

@everywhere include("KuruczKorgBridge.jl")
@everywhere using .KuruczKorgBridge

# Distribute work over workers
results = pmap(params -> generate_spectrum_for_params(emulator, params, wavelengths),
               param_grid)
```

---

## References

### Kurucz Models

- Li, J. (2025). *Kurucz-a1: Kurucz Stellar Atmosphere Emulator*. GitHub: [jiadonglee/kurucz1](https://github.com/jiadonglee/kurucz1)
- Kurucz, R. L. (1993). *ATLAS9 Stellar Atmosphere Programs and 2 km/s Grid*. Kurucz CD-ROM No. 13.

### Korg.jl

- Wheeler, A., Kochukhov, O., et al. (2024). *Korg.jl: Efficient 1D LTE Spectral Synthesis in Julia*.
- Documentation: [Korg.jl Docs](https://ajwheeler.github.io/Korg.jl/stable/)

### JWST NIRSpec

- NIRSpec Instrument Handbook: [STScI Documentation](https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph)
- Wavelength ranges, resolution, and sensitivity curves

### M31 Stellar Spectroscopy

- Relevantpapers for your specific science case (add as appropriate)

---

## Contributing & Support

For issues specific to this integration:
1. Check this README first
2. Verify installation and configuration
3. Test with example scripts
4. Check error messages in `failed_params.txt`

For issues with:
- **Kurucz emulator**: [jiadonglee/kurucz1 issues](https://github.com/jiadonglee/kurucz1/issues)
- **Korg.jl**: [ajwheeler/Korg.jl issues](https://github.com/ajwheeler/Korg.jl/issues)

---

## Citation

If you use this pipeline in your research, please cite:

1. **Kurucz Emulator:**
   ```
   Li, J. (2025). Physics-Informed Neural Emulation of Kurucz Stellar
   Atmospheric Models with Optical Depth Integration.
   GitHub: github.com/jiadonglee/kurucz1
   ```

2. **Korg.jl:**
   ```
   Wheeler, A., et al. (2024). Korg.jl: Efficient 1D LTE Spectral
   Synthesis in Julia. [Add paper citation when published]
   ```

3. **Original Kurucz Models:**
   ```
   Kurucz, R. L. (1993). ATLAS9 Stellar Atmosphere Programs and 2 km/s grid.
   Kurucz CD-ROM No. 13. Cambridge, Mass.: Smithsonian Astrophysical Observatory.
   ```

---

## License

This integration code is provided under the same license as Korg.jl (MIT License).
See the main Korg.jl LICENSE file for details.

---

**Last updated:** 2025-10-24
**Version:** 1.0
**Contact:** See Korg.jl maintainers
