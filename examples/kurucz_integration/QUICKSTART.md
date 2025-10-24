# Quick Start Guide: Kurucz-Korg Pipeline

**Goal:** Generate synthetic stellar spectra for JWST M31 analysis in under 10 minutes.

---

## 1. Prerequisites Checklist

- [ ] Julia 1.6+ installed (`julia --version`)
- [ ] Python 3.8+ installed (`python --version`)
- [ ] Git installed

---

## 2. Installation (5 minutes)

### Step 1: Install Kurucz Emulator

```bash
# Clone and install
git clone https://github.com/jiadonglee/kurucz1.git
cd kurucz1
pip install -e .

# Verify
python -c "import kuruczone; print('✓ Kurucz installed')"
```

### Step 2: Setup Julia Environment

```bash
cd /path/to/Korg.jl/examples/kurucz_integration
julia
```

In Julia REPL:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()  # Installs all dependencies

# Configure PyCall
using PyCall
# If needed: ENV["PYTHON"] = "/path/to/python"
# Pkg.build("PyCall")

# Test
pyimport("kuruczone")
using Korg
```

### Step 3: Download Model Weights

Download `a_one_weights.pt` and set:
```bash
export KURUCZ_MODEL_PATH="/path/to/a_one_weights.pt"
```

---

## 3. Generate Your First Spectrum (2 minutes)

```julia
# In Julia REPL
include("example_single_spectrum.jl")
```

**Output:**
- `spectrum_5777K_logg4.44_FeH0.0.dat` - Spectrum data file
- `spectrum_plot.png` - Visualization (if Plots.jl installed)

**What just happened?**
1. Loaded Kurucz emulator
2. Generated solar atmosphere (Teff=5777K)
3. Converted to Korg format
4. Synthesized JWST NIRSpec spectrum (1.5-2.8 μm)

---

## 4. Generate a Test Training Grid (3 minutes)

```julia
include("generate_jwst_training_grid.jl")
```

**Output:**
- `jwst_m31_training_spectra.h5` - HDF5 file with 5 test spectra

**Grid contains:**
- Sun (Teff=5777K, logg=4.44, [Fe/H]=0.0)
- Metal-poor dwarf (5000K, 4.5, -0.5)
- Metal-poor giant (4500K, 2.5, -1.0)
- Hot dwarf (6000K, 4.0, 0.0)
- Cool giant (4000K, 1.5, -0.5)

---

## 5. Inspect Results

### In Julia:
```julia
using HDF5

h5open("jwst_m31_training_spectra.h5", "r") do f
    println("Number of models: ", read(f["metadata/n_models"]))
    println("Wavelength range: ", extrema(read(f["wavelengths"])))
    println("Parameters:")
    println("  Teff: ", read(f["parameters/Teff"]))
    println("  logg: ", read(f["parameters/logg"]))
end
```

### In Python:
```python
import h5py
import matplotlib.pyplot as plt

with h5py.File('jwst_m31_training_spectra.h5', 'r') as f:
    wl = f['wavelengths'][:]
    spectra = f['spectra'][:]
    Teff = f['parameters/Teff'][:]

    # Plot first spectrum
    plt.plot(wl, spectra[0])
    plt.xlabel('Wavelength (Å)')
    plt.ylabel('Flux')
    plt.title(f'Teff = {Teff[0]:.0f} K')
    plt.show()
```

---

## 6. Generate Production Grid

**Edit `generate_jwst_training_grid.jl`:**

```julia
# Change line:
use_test_grid = false  # Was: true

# Optionally adjust parameter ranges:
param_grid = get_stellar_parameter_grid(
    Teff_min=3500.0, Teff_max=7000.0, Teff_step=250.0,
    logg_min=1.0, logg_max=5.0, logg_step=0.5,
    FeH_min=-2.5, FeH_max=0.5, FeH_step=0.5,
    alphaFe_min=0.0, alphaFe_max=0.4, alphaFe_step=0.2
)
```

**Run:**
```julia
include("generate_jwst_training_grid.jl")
```

**⚠️ Warning:** This will generate ~2500 spectra and take 10-20 hours on a single core!

**For faster results:**
- Increase step sizes (e.g., `Teff_step=500.0`)
- Reduce wavelength range
- Implement parallelization (see README)

---

## 7. Common Issues & Quick Fixes

| Issue | Quick Fix |
|-------|-----------|
| `ERROR: PyError` | Check: `python -c "import kuruczone"` |
| `KURUCZ_MODEL_PATH not found` | Set: `ENV["KURUCZ_MODEL_PATH"] = "/path/to/model.pt"` |
| Synthesis fails | Check stellar parameters are in valid range |
| Very slow | Increase wavelength step size (e.g., 0.5 → 1.0 Å) |
| Out of memory | Reduce grid size or wavelength range |

---

## 8. Next Steps

✅ **You're ready to:**
1. Customize parameter grid for your science case
2. Adjust wavelength ranges for your JWST program
3. Generate production training dataset
4. Train ML models for parameter inference
5. Apply to JWST M31 observations

📖 **For detailed information, see:**
- `README.md` - Full documentation with limitations
- `KuruczKorgBridge.jl` - API documentation
- `generate_jwst_training_grid.jl` - Customization options

---

## Questions?

- **Kurucz emulator issues:** https://github.com/jiadonglee/kurucz1/issues
- **Korg.jl issues:** https://github.com/ajwheeler/Korg.jl/issues
- **This pipeline:** Check README.md troubleshooting section

---

**Estimated time to first spectrum:** < 10 minutes
**Estimated time to production grid:** Hours to days (depending on size)

Good luck with your JWST M31 analysis! 🚀🔭
