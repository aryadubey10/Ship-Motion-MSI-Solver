# 🚢 Ship Motion & Motion Sickness Incidence (MSI) Solver

A **six-degree-of-freedom (6-DOF) frequency-domain seakeeping solver** that predicts Motion Sickness Incidence (MSI) aboard ships using the **O'Hanlon & McCauley model**. Built in Python with NumPy/SciPy.

![Python](https://img.shields.io/badge/Python-3.10+-blue?logo=python&logoColor=white)
![NumPy](https://img.shields.io/badge/NumPy-2.1-013243?logo=numpy)
![SciPy](https://img.shields.io/badge/SciPy-1.14-8CAAE6?logo=scipy)
![License](https://img.shields.io/badge/License-MIT-green)

---

## 🎯 What It Does

Given hydrodynamic input data (added mass, radiation damping, wave excitation forces), this solver:

1. **Solves the coupled 6-DOF equations of motion** in the frequency domain via batched complex linear algebra
2. **Computes JONSWAP wave spectra** (IACS Rec. 34 normalised) for realistic North Atlantic sea states
3. **Calculates cabin vertical acceleration spectra** including encounter-frequency transformation and short-wave cosine tapering
4. **Predicts Motion Sickness Incidence** (% of passengers vomiting within 2 hours) using the O'Hanlon & McCauley statistical model
5. **Determines roll natural frequency** via fixed-point iteration with PCHIP-interpolated added mass
6. **Performs sensitivity studies** on cabin position, damping ratio, wave heading, and ship speed

---

## 📐 Physics & Models

| Component | Method |
|:--|:--|
| Wave Spectrum | IACS Rec. 34 JONSWAP (Hₛ = 4.5 m, T₀₋₁ = 11.5 s, γ = 1.5) |
| Equations of Motion | 6-DOF coupled: `[-ωₑ²(M + A) + iωₑB + C] ξ̂ = F̂` |
| Roll Damping | Wave-making (from radiation) + viscous (ζᵥ = 0.05 of critical) |
| Roll Natural Frequency | Fixed-point iteration with Brent root-finding |
| Interpolation | PCHIP (roll region) / Cubic Spline (fine grid) |
| Short-Wave Filter | Cosine taper beyond `√(8πg/L)` for head & bow-quartering seas |
| MSI Model | O'Hanlon & McCauley: `MSI = 100 · Φ[(log₁₀(ā/g) − μ) / 0.4]` |
| Spectral Moments | Trapezoidal integration on wave-frequency axis |

---

## 🗂️ Project Structure

```
Ship-Motion-MSI-Solver/
├── src/
│   ├── msi.py                    # Main solver & CLI entry point
│   ├── import_data.py            # HydRA data ingestion
│   ├── build_report.py           # Report/figure generation
│   └── msi_project/              # Core library modules
│       ├── spectrum.py           # JONSWAP spectrum & spectral moments
│       ├── motions.py            # 6-DOF solve, encounter freq, roll damping
│       ├── msi.py                # O'Hanlon & McCauley MSI model
│       ├── cabin.py              # Cabin position & vertical acceleration
│       ├── clean.py              # Irregular-frequency spike removal
│       ├── io_hydra.py           # HydRA JSON I/O
│       ├── cases.py              # Case matrix generation
│       ├── plots.py              # Visualization utilities
│       └── constants.py          # Physical & project constants
├── results/                      # Computed output CSVs
│   ├── msi_cases.csv             # 30-case MSI results (2 hulls × 3 headings × 5 speeds)
│   ├── roll_damping.csv          # Roll natural frequency & damping coefficients
│   ├── sensitivity.csv           # Parameter sensitivity studies
│   └── spikes.csv                # Irregular-frequency spike log
├── report/
│   ├── report.pdf                # Technical report
│   └── figures/                  # Generated plots
├── requirements.txt
├── pytest.ini
└── LICENSE
```

---

## 🚀 Quick Start

```bash
# Clone the repository
git clone https://github.com/aryadubey10/Ship-Motion-MSI-Solver.git
cd Ship-Motion-MSI-Solver

# Set up the environment
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Run the solver (default variant: KCS_L0p950)
python src/msi.py KCS_L0p950

# Custom integration step for convergence studies
python src/msi.py KCS_L0p950 --dw 0.0005
```

---

## 📊 Output

The solver produces:

- **30-case MSI table** (`msi_cases.csv`): σₐ, ā, ω̄ₑ, μ_MSI, MSI% for every hull/heading/speed combination
- **Roll damping parameters** (`roll_damping.csv`): ωₙ, ζ_wave, B_v,44
- **Sensitivity studies** (`sensitivity.csv`): cabin position, viscous damping ratio effects
- **Publication-quality figures**: wave spectra, roll RAOs, acceleration spectra, MSI contour maps

### Sample Figures

| MSI vs Speed | Acceleration Spectrum | MSI Contour Map |
|:--:|:--:|:--:|
| MSI at three wave headings | Cabin acceleration in encounter-freq domain | All 30 cases on equal-MSI contours |

---

## 🔧 Technical Highlights

- **Batched linear solve**: Vectorised `np.linalg.solve` over all frequencies — no Python loops in the hot path
- **PCHIP vs Cubic Spline**: PCHIP used for roll-region interpolation (0.1–0.2 rad/s) where B₄₄/A₄₄ vary steeply; cubic spline for smooth fine-grid hydrodynamics
- **Brent root-finding**: Robust roll natural-frequency computation with bracket detection
- **Encounter-frequency Jacobian**: Proper `dω/dωₑ` transformation for spectral density on encounter axis
- **Validation pipeline**: Zero-viscous RAOs verified against HydRA reference solutions (< 0.1% relative error)

---

## 📚 References

- O'Hanlon, J.F. and McCauley, M.E. (1974). *Motion Sickness Incidence as a Function of the Frequency and Acceleration of Vertical Sinusoidal Motion*. Aerospace Medicine.
- IACS Recommendation No. 34 (2001). *Standard Wave Data*.
- Journée, J.M.J. and Massie, W.W. (2001). *Offshore Hydromechanics*. Delft University of Technology.

---

## 📄 License

MIT License — see [LICENSE](LICENSE) for details.

---

*Built as part of naval architecture research on passenger comfort and hull-form optimisation.*
