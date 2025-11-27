# Kurucz-Korg Pipeline: Complete Workflow Summary

## Overview

This document provides a concise summary of the complete workflow for generating JWST training spectra using Kurucz atmospheres and Korg synthesis.

---

## Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                    STELLAR PARAMETERS                           │
│         Teff (K), logg, [Fe/H], [α/Fe]                         │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│              KURUCZ NEURAL EMULATOR (Python)                    │
│  Input: [Teff, logg, [Fe/H], [α/Fe], τ_grid]                  │
│  Output: T(τ), P(τ), ρ(τ), n_e(τ), κ(τ), ...                 │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│              CONVERSION LAYER (Julia)                           │
│  • Calculate geometric heights: z(τ)                            │
│  • Calculate number density: n = P/(k_B*T)                      │
│  • Create Korg atmosphere layers                                │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│              KORG SPECTRAL SYNTHESIS (Julia)                    │
│  Input: Atmosphere, abundances, λ range, linelist              │
│  Output: F_λ (emergent spectrum)                               │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│              JWST WAVELENGTH SAMPLING                           │
│  NIRSpec: 1.0-5.2 μm (10000-52000 Å)                          │
└────────────────────────┬────────────────────────────────────────┘
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│              TRAINING DATASET (HDF5)                            │
│  Parameters × Spectra → ML Training                             │
└─────────────────────────────────────────────────────────────────┘
```

---

## Step-by-Step Workflow

### Phase 1: Setup

1. **Install Kurucz emulator** (Python)
   ```bash
   git clone https://github.com/jiadonglee/kurucz1.git
   cd kurucz1
   pip install -e .
   ```

2. **Download model weights**
   - Get `a_one_weights.pt`
   - Set `KURUCZ_MODEL_PATH` environment variable

3. **Configure Julia environment**
   ```julia
   using Pkg
   Pkg.add(["Korg", "PyCall", "HDF5", "ProgressMeter"])
   ```

### Phase 2: Single Spectrum Generation

**File:** `example_single_spectrum.jl`

```julia
# Load emulator
emulator = load_kurucz_emulator(model_path)

# Generate atmosphere
kurucz_atm = get_kurucz_atmosphere(emulator, Teff, logg, Fe_H, alpha_Fe)

# Convert to Korg format
korg_atm = kurucz_to_korg(kurucz_atm)

# Synthesize spectrum
A_X = Korg.format_A_X(Fe_H, alpha_H=alpha_Fe)
result = synthesize_spectrum(korg_atm, A_X, wavelengths; vmic=1.0)
```

**Or use convenience function:**
```julia
result = kurucz_synthesize(emulator, Teff, logg, Fe_H, alpha_Fe, wavelengths)
```

### Phase 3: Grid Generation

**File:** `generate_jwst_training_grid.jl`

```julia
# Define parameter grid
param_grid = get_stellar_parameter_grid(
    Teff_min=3500.0, Teff_max=7000.0, Teff_step=250.0,
    logg_min=1.0, logg_max=5.0, logg_step=0.5,
    FeH_min=-2.5, FeH_max=0.5, FeH_step=0.5,
    alphaFe_min=0.0, alphaFe_max=0.4, alphaFe_step=0.2
)

# Choose JWST wavelength configuration
wavelengths = JWST_WAVELENGTH_CONFIGS[:nirspec_g235m]

# Generate all spectra
results = generate_grid(emulator, param_grid, wavelengths;
                       output_file="jwst_training.h5")
```

### Phase 4: Data Usage (Python/Julia)

**Load training data:**
```python
import h5py
with h5py.File('jwst_training.h5', 'r') as f:
    X = f['spectra'][:]           # Features
    y_Teff = f['parameters/Teff'][:]    # Labels
    y_logg = f['parameters/logg'][:]
    y_FeH = f['parameters/FeH'][:]
```

**Train ML model:**
```python
from sklearn.ensemble import RandomForestRegressor
# or use neural networks, etc.

model = RandomForestRegressor()
model.fit(X, np.column_stack([y_Teff, y_logg, y_FeH]))
```

**Apply to JWST observations:**
```python
# Load JWST M31 spectrum
jwst_spectrum = load_jwst_spectrum("m31_star.fits")

# Predict parameters
predicted_params = model.predict(jwst_spectrum)
Teff, logg, FeH = predicted_params[0]
```

---

## File Structure

```
examples/kurucz_integration/
├── README.md                              # Full documentation
├── QUICKSTART.md                          # 10-minute getting started
├── WORKFLOW_SUMMARY.md                    # This file
├── requirements.txt                       # Python dependencies
├── Project.toml                           # Julia dependencies
├── KuruczKorgBridge.jl                   # Core integration module
├── generate_jwst_training_grid.jl        # Grid generation script
└── example_single_spectrum.jl            # Single spectrum example
```

---

## Key Functions Reference

### Loading & Setup

| Function | Purpose | Usage |
|----------|---------|-------|
| `load_kurucz_emulator(path)` | Load neural net model | `emulator = load_kurucz_emulator("model.pt")` |

### Atmosphere Generation

| Function | Purpose | Returns |
|----------|---------|---------|
| `get_kurucz_atmosphere(...)` | Generate atm structure | `KuruczAtmosphere` struct |
| `kurucz_to_korg(kurucz_atm)` | Convert to Korg format | `PlanarAtmosphere` |

### Synthesis

| Function | Purpose | Returns |
|----------|---------|---------|
| `synthesize_spectrum(...)` | Korg synthesis | `(wavelengths, flux, continuum)` |
| `kurucz_synthesize(...)` | End-to-end | Spectrum results |

### Grid Generation

| Function | Purpose | Returns |
|----------|---------|---------|
| `get_stellar_parameter_grid(...)` | Create param grid | Vector of NamedTuples |
| `generate_grid(...)` | Batch synthesis | Saves HDF5 file |

---

## Data Flow

### Input Data
- **Stellar parameters:** 4 values per star
- **Optical depth grid:** 80 points (default)
- **Wavelength grid:** User-defined (JWST ranges)

### Intermediate Data (per star)
- **Kurucz atmosphere:** 7 arrays × 80 depths = 560 floats
- **Korg atmosphere:** Converted structure + derived quantities

### Output Data
- **Single spectrum:** ~26,000 wavelength points (1.5-2.8 μm @ 0.5 Å)
- **Grid of 2500 spectra:** ~200 MB HDF5 file

---

## Performance Benchmarks

*Approximate timings on typical workstation (single core):*

| Task | Time | Notes |
|------|------|-------|
| Load Kurucz emulator | ~5 sec | One-time per session |
| Generate Kurucz atmosphere | ~0.1 sec | 80 depth points |
| Convert to Korg | ~0.01 sec | Fast |
| Load linelist | ~5 sec | One-time per session |
| Synthesize spectrum (1.5-2.8 μm @ 0.5 Å) | ~15 sec | Depends on line density |
| **Total per star** | ~20 sec | After initial loading |
| Grid of 2500 stars | ~14 hours | Can parallelize |

**Speedup strategies:**
1. Coarser wavelength sampling (0.5 → 1.0 Å): 2× faster
2. Narrower wavelength range: Linear speedup
3. Parallelization (8 cores): ~7× faster
4. Reduce grid size: Linear speedup

---

## Parameter Ranges

### Recommended Grid (based on M31 stellar populations)

| Parameter | Min | Max | Step | N points |
|-----------|-----|-----|------|----------|
| Teff (K) | 3500 | 7000 | 250 | 15 |
| logg | 1.0 | 5.0 | 0.5 | 9 |
| [Fe/H] | -2.5 | 0.5 | 0.5 | 7 |
| [α/Fe] | 0.0 | 0.4 | 0.2 | 3 |
| **Total** | | | | **2835** |

### Physical Interpretation

- **Teff 3500-5000K:** RGB, AGB stars
- **Teff 5000-6500K:** Main sequence turn-off, subgiants
- **Teff 6500-7000K:** Hot main sequence
- **logg 1.0-2.5:** Giants, supergiants
- **logg 3.0-4.0:** Subgiants
- **logg 4.5-5.0:** Dwarfs
- **[Fe/H] < -1.5:** Halo population
- **[Fe/H] -1.0 to -0.5:** Thick disk
- **[Fe/H] > -0.5:** Thin disk
- **[α/Fe] 0.2-0.4:** Old, rapid formation

---

## JWST Wavelength Configurations

| Config | λ Range (μm) | λ Range (Å) | Use Case |
|--------|-------------|-------------|----------|
| NIRSpec Prism | 0.6-5.3 | 6000-53000 | Full coverage, low-res |
| NIRSpec G140M | 0.97-1.84 | 9700-18400 | H-band, medium-res |
| NIRSpec G235M | 1.7-3.07 | 17000-30700 | K-band, medium-res |
| NIRSpec G395M | 2.9-5.18 | 29000-51800 | L/M-band, medium-res |
| **M31 Survey** | **1.5-2.8** | **15000-28000** | **Stellar spectroscopy** |

---

## Critical Reminders

⚠️ **Before running production grid:**

1. ✅ Verify Kurucz model training range
2. ✅ Test with small grid first (`use_test_grid = true`)
3. ✅ Check atmospheric structure for edge cases
4. ✅ Validate against known stars (Sun, Arcturus, etc.)
5. ✅ Ensure sufficient disk space (~200 MB per 2500 spectra)
6. ✅ Consider parallelization for large grids

⚠️ **Physical limitations:**

- 1D plane-parallel atmospheres
- LTE assumption (no NLTE)
- No rotation, macroturbulence
- Line list incompleteness (cool stars, IR)

⚠️ **For M31 science:**

- Apply extinction corrections separately
- Account for different abundance patterns
- Consider stellar multiplicity in crowded fields
- Include radial velocity shifts

---

## Quick Command Reference

```julia
# Single spectrum (solar)
result = kurucz_synthesize(emulator, 5777, 4.44, 0.0, 0.0, [(15000, 28000, 0.5)])

# Test grid (5 stars)
include("generate_jwst_training_grid.jl")  # with use_test_grid=true

# Production grid
# Edit use_test_grid=false, then:
include("generate_jwst_training_grid.jl")

# Load results
using HDF5
h5open("jwst_training.h5", "r") do f
    Teff = read(f["parameters/Teff"])
    spectra = read(f["spectra"])
end
```

---

## Troubleshooting Checklist

- [ ] Python kuruczone installed: `python -c "import kuruczone"`
- [ ] PyCall configured: `using PyCall; pyimport("torch")`
- [ ] Model path set: `ENV["KURUCZ_MODEL_PATH"]`
- [ ] Korg installed: `using Korg`
- [ ] HDF5 installed: `using HDF5`
- [ ] Parameters in valid range (see README)
- [ ] Sufficient disk space
- [ ] Tested with single spectrum first

---

## Next Steps

1. **Validation:** Compare synthetic spectra to observations
2. **Training:** Use HDF5 output for ML model training
3. **Application:** Apply trained model to JWST M31 data
4. **Iteration:** Refine grid based on initial results
5. **Science:** Stellar populations, chemical evolution, etc.

---

**For full details, see README.md**
**For quick start, see QUICKSTART.md**
**For API reference, see KuruczKorgBridge.jl**

---

Generated: 2025-10-24
Version: 1.0
Branch: `claude/kurucz-stellar-spectra-generation-011CURCdKUaabLyzvzNdrehT`
