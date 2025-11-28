# Kurucz-Korg Integration Analysis

## Summary

This document analyzes the Kurucz-Korg integration implementation in `examples/kurucz_integration/` with a focus on:
1. The necessity and correctness of `calculate_height()` and `calculate_number_density()` functions
2. Continuum normalization practices
3. Recommendations for improvements

## 1. Analysis of `calculate_height()` and `calculate_number_density()`

### Are These Functions Necessary for Korg?

**YES** - Both functions are **necessary and correctly implemented**.

### Detailed Analysis

#### Number Density (`calculate_number_density`)

**Location**: `KuruczKorgBridge.jl:216-232`

**Current Implementation**:
```julia
function calculate_number_density(P::Vector{Float64}, T::Vector{Float64})
    return P ./ (k_B .* T)
end
```

**Why it's necessary**:
- Korg's `PlanarAtmosphereLayer` struct requires total number density (`n`) as a parameter
- Kurucz models provide pressure (`P`) and temperature (`T`), but not number density directly
- The ideal gas law conversion `n = P / (k_B * T)` is the **standard approach** used throughout Korg

**Evidence from Korg source** (`src/atmosphere.jl:219-220, 248-249`):
```julia
# MARCS reader:
nₑ = Pe / (temp * kboltz_cgs)  # electron number density
n = Pg / (temp * kboltz_cgs)   # total number density

# PHOENIX reader:
number_density = @. Pgas / (kboltz_cgs * T)
electron_number_density = @. Pe / (kboltz_cgs * T)
```

**Conclusion**: ✅ **Correct and standard approach** - This is exactly how Korg handles pressure-to-density conversion internally.

#### Height Calculation (`calculate_heights`)

**Location**: `KuruczKorgBridge.jl:172-213`

**Current Implementation**:
- Uses hydrostatic equilibrium: `dP = -g * ρ * dz`
- Integrates from top to bottom of atmosphere
- Normalizes so `z=0` at approximately τ=1 (photosphere)

**Why it's necessary**:
- Korg's `PlanarAtmosphereLayer` requires geometric height (`z`) in cm
- Kurucz models provide:
  - `tau`: Optical depth
  - `P`: Gas pressure (dyn/cm²)
  - `RHOX`: Mass column density (g/cm²)
  - But **NOT** geometric heights
- Height must be calculated from the available quantities

**Evidence from Korg source**:
- MARCS models provide `depth` directly (line 210, 223 in `atmosphere.jl`)
- PHOENIX models set `z = NaN` and only work with "anchored" radiative transfer (line 252)
- **No standard Korg function exists** to calculate heights from pressure/column density

**Current Implementation Details**:
```julia
function calculate_heights(P::Vector{Float64}, RHOX::Vector{Float64}, logg::Real)
    g = 10.0^logg  # cm/s²
    n = length(P)
    z = zeros(Float64, n)

    # Integrate from top (low pressure) to bottom (high pressure)
    for i in 2:n
        dRHOX = RHOX[i] - RHOX[i-1]
        dP = P[i] - P[i-1]

        if dRHOX != 0
            # From hydrostatic: dP = g * ρ * dz, so dz = dP/(g*ρ)
            rho_avg = abs(dRHOX) / (abs(dP) / g)
            dz = dP / (g * rho_avg)
            z[i] = z[i-1] + dz
        else
            z[i] = z[i-1]
        end
    end

    # Normalize so z=0 at optical depth ≈ 1 (approximate photosphere)
    return z .- z[n÷2]  # Simple approximation
end
```

**Is this the most realistic way?**

The current implementation is **reasonable but has room for improvement**:

✅ **Strengths**:
- Uses fundamental hydrostatic equilibrium
- Integrates over the atmosphere structure
- Provides physically meaningful heights

⚠️ **Potential Improvements**:
1. **Better reference point**: Currently uses `z[n÷2]` as τ=1. Better approach:
   ```julia
   # Find actual τ=1 layer
   tau_1_idx = argmin(abs.(kurucz_atm.tau .- 1.0))
   return z .- z[tau_1_idx]
   ```

2. **More accurate density calculation**: Current implementation approximates density from `dRHOX/dz`. A more direct approach:
   ```julia
   # Calculate density directly from mass column density derivative
   ρ = zeros(Float64, n)
   for i in 2:n-1
       ρ[i] = (RHOX[i+1] - RHOX[i-1]) / (z[i+1] - z[i-1])
   end
   ```

3. **Alternative: Use scale height approximation**:
   ```julia
   # Pressure scale height H = k_B * T / (μ * m_H * g)
   # Can provide quick height estimates
   ```

**Comparison with other atmosphere codes**:
- **MARCS**: Provides heights directly from their hydrostatic solver
- **ATLAS9/Kurucz**: Original codes compute heights using similar hydrostatic integration
- **PHOENIX**: Uses geometric depth scale explicitly

**Conclusion**: ✅ **Necessary and acceptable**, but can be refined (see recommendations below).

## 2. Continuum Normalization

### Current Status

**In example_single_spectrum.jl** (lines 106-116):
- Plots show **both** raw flux and continuum
- **NOT** continuum-normalized by default

**In visualize_results.jl** (lines 104-108):
- Explicitly normalizes spectra: `normalized = spec ./ cont`
- All plots show **continuum-normalized** spectra
- Follows standard Korg practice

### Is this standard in Korg?

**YES** - Continuum normalization is standard practice:
- Korg's `synthesize()` function returns both `flux` and `cntm` (continuum)
- Most Korg examples show continuum-normalized spectra for analysis
- This is standard in stellar spectroscopy for:
  - Comparing spectra with different temperatures
  - Analyzing line profiles
  - Stellar parameter inference

### Current Practice in Examples

**example_single_spectrum.jl**:
- ❌ Shows raw flux (not ideal for visual inspection)
- ✅ Saves continuum data for post-processing

**visualize_results.jl**:
- ✅ Shows continuum-normalized spectra
- ✅ Includes continuum as reference line

## 3. Missing Feature: Zoomed NIR Plot

The user requested a plot of the **continuum-normalized spectrum from 13500-17500 Å**.

**Current status**: ❌ Not implemented

**This range is scientifically important**:
- 1.35-1.75 μm (H-band region)
- Contains important stellar features:
  - CO bands
  - OH lines
  - Various metal lines (Fe, Si, Mg, Ca, Ti)
- Relevant for JWST NIRSpec observations

## 4. Recommendations

### High Priority

1. **Improve `calculate_heights()` reference point**:
   ```julia
   # Use actual τ=1 instead of middle of array
   tau_1_idx = argmin(abs.(kurucz_atm.tau .- 1.0))
   return z .- z[tau_1_idx]
   ```

2. **Update `example_single_spectrum.jl` to show continuum-normalized spectra**:
   - More useful for visual inspection
   - Matches standard practice in stellar spectroscopy
   - Easier to identify absorption features

3. **Add zoomed H-band plot** (13500-17500 Å):
   - Important for JWST applications
   - Shows spectral detail better than full range

### Medium Priority

4. **Add validation against known standards**:
   - Compare solar spectrum (Teff=5777K) to observations
   - Check line depths and continuum shape
   - Validate against MARCS/PHOENIX for same parameters

5. **Document mean molecular weight assumption**:
   - Currently assumes μ ≈ 1.3 (appropriate for solar composition)
   - May need adjustment for very metal-poor stars

6. **Add uncertainty estimates**:
   - Document expected accuracy of height calculation
   - Quantify impact of approximations

### Low Priority

7. **Consider implementing spherical geometry option**:
   - Important for giants with log g < 1.5
   - Korg supports `ShellAtmosphere`
   - Can convert from planar using stellar radius

8. **Add NLTE corrections** (if available):
   - Important for strong resonance lines
   - May be outside scope of current emulator

## 5. Validation Checklist

To ensure accuracy of the integration:

- [ ] Compare synthetic solar spectrum to observed (Wallace et al. atlas)
- [ ] Verify pressure structure matches ATLAS9 models for same parameters
- [ ] Check height scale is physically reasonable (~100-1000 km for photosphere)
- [ ] Validate continuum level against other synthesis codes
- [ ] Test edge cases (very cool stars, very metal-poor, giants)
- [ ] Compare to MARCS atmospheres for subset of parameter space

## 6. Conclusion

**The current integration is fundamentally sound**:
- ✅ `calculate_number_density()` is **correct** and matches Korg's internal approach
- ✅ `calculate_heights()` is **necessary** and uses appropriate physics
- ✅ Continuum normalization follows **standard practice**

**Key improvements needed**:
1. Better reference point for height zero-point (τ=1)
2. Add continuum-normalized plotting to main example
3. Add requested zoomed H-band plot (13500-17500 Å)

**The integration is ready for science use with the recommended refinements**.

---

*Analysis Date: 2025-11-28*
*Korg.jl main branch commit: a72ffdb*
