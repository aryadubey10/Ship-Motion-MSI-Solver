# Report — OE3035 Project: Motion Sickness Incidence on Board

**Group 10 · KCS_L0p950 vs KCS_L1p000**
Course: OE3035 Motion of Ships and Floating Systems · Department of Ocean Engineering, IIT Madras
Instructor: Prof. Abhilash Somayajula

---

## 📄 Submitted Document

| File | Size | Description |
|:--|:--|:--|
| [`report.pdf`](report.pdf) | 802 KB | Full project report (Deliverables 1–5) |

The PDF was generated from `src/build_report.py` using the numbers written to `results/` by `src/msi.py`. Every number in the report traces directly to a computed result — nothing is hand-typed.

---

## 📊 Figures in This Folder (`figures/`)

| File | Deliverable | Description |
|:--|:--|:--|
| [`figures/roll_rao.png`](figures/roll_rao.png) | Fig. 1 | Roll RAO \|ξ₄\| (rad/m) vs ω_e in beam seas (β = 90°) at Fn_des, both hulls, with and without viscous roll damping B_v,44. ω_n marked. |
| [`figures/msi_vs_speed.png`](figures/msi_vs_speed.png) | Fig. 2 | MSI (%) vs Fn/Fn_des, one line per heading (90°/135°/180°), both hulls overlaid. |
| [`figures/accel_spectrum_90.png`](figures/accel_spectrum_90.png) | Fig. 3a | Vertical acceleration spectrum S_aa(ω_e) at Fn_des — beam seas (β = 90°). |
| [`figures/accel_spectrum_135.png`](figures/accel_spectrum_135.png) | Fig. 3b | Vertical acceleration spectrum S_aa(ω_e) at Fn_des — bow-quartering seas (β = 135°). |
| [`figures/accel_spectrum_180.png`](figures/accel_spectrum_180.png) | Fig. 3c | Vertical acceleration spectrum S_aa(ω_e) at Fn_des — head seas (β = 180°). |
| [`figures/msi_contours.png`](figures/msi_contours.png) | Fig. 4 | All 30 cases on the (ω̄_e, ā) plane with MSI contours at 2, 5, 10, 15, 20, 30%. Marker shape = heading, fill = hull, size = speed. |
| [`figures/sea_spectrum.png`](figures/sea_spectrum.png) | — | JONSWAP wave spectrum S(ω) for the IACS Rec. 34 North Atlantic design condition. |

Additional deliverable figures (PDF + 300-dpi PNG) with improved styling are in [`Extra Research/figures/`](../Extra%20Research/figures/).

---

## 📋 Deliverables Checklist

### Deliverable 1 — Derivation of the Cabin Acceleration Spectrum
*(Report §4: Method)*

The derivation proceeds step by step:

1. **Regular wave → 6-DOF EoM** (frequency domain, Ch. 8):

   Z(ω_e) ξ̂ = F̂_FK + F̂_S,   Z = −ω_e²(M + A) + iω_e(B + B_v) + C

2. **Encounter frequency** (deep water, Ch. 7):

   ω_e = ω − (ω² U / g) cos β

3. **Linearised vertical displacement at cabin point P** (Ch. 4, δz = φy − θx):

   ŝ₃(ω) = ξ̂₃ + y_P ξ̂₄ − x_P ξ̂₅

4. **Vertical acceleration RAO** (d/dt → iω_e, Ch. 3):

   Ĥ_a(ω) = −ω_e² ŝ₃(ω)

5. **Response spectrum integrated over wave frequency** (Ch. 9):

   m_n^a = ∫ |ω_e(ω)|ⁿ |Ĥ_a(ω)|² S(ω) dω,   n = 0, 1

6. **Statistical measures**:

   σ_a = √(m₀^a),   ω̄_e = m₁^a / m₀^a,   ā = √(2/π) σ_a

**Assumptions stated:** linear potential flow, small-amplitude motions, stationary ergodic Gaussian sea, long-crested waves, deep water, point cabin, linearised roll damping, no slamming/green water/parametric roll.

---

### Deliverable 2 — Roll Damping

*(Report §6, Table 1; computed in `results/roll_damping.csv`)*

| Hull | ω_n (rad/s) | Roll period T_n (s) | ζ_wave | B_v,44 (N·m·s) |
|:--|:--|:--|:--|:--|
| KCS_L1p000 | **0.164596** | 38.17 | 3.543 × 10⁻⁵ | **1.8192 × 10⁸** |
| KCS_L0p950 | **0.162608** | 38.64 | 3.351 × 10⁻⁵ | **1.8895 × 10⁸** |

ω_n found by fixed-point iteration using A₄₄ at β = 90°, Fn_des. The tiny ζ_wave (≈ 3.5 × 10⁻⁵) confirms that potential flow provides almost no roll damping — essentially all damping comes from the added viscous term ζ_v = 0.05.

**Irregular-frequency spikes** are documented in `results/spikes.csv` (27 replacement events across both hulls, mostly at high encounter frequencies in 135° and 180° headings). The detection rule and replaced values are listed in the report's verification section.

---

### Deliverable 3 — 30-Case Results Table

*(Report §6, Table 2; full data in `results/msi_cases.csv`)*

**KCS_L1p000 (baseline, Lpp = 230.0 m)**

| β (°) | Fn/Fn_des | U (m/s) | σ_a (m/s²) | ā (m/s²) | ω̄_e (rad/s) | MSI (%) |
|:--|:--|:--|:--|:--|:--|:--|
| 90  | 0.90 | 11.1132 | 0.4337 | 0.3460 | 0.6716 | **3.947** |
| 90  | 0.95 | 11.7306 | 0.4349 | 0.3470 | 0.6716 | **3.972** |
| 90  | 1.00 | 12.3480 | 0.4360 | 0.3479 | 0.6716 | **3.997** |
| 90  | 1.05 | 12.9654 | 0.4372 | 0.3489 | 0.6716 | **4.023** |
| 90  | 1.10 | 13.5828 | 0.4384 | 0.3498 | 0.6717 | **4.049** |
| 135 | 0.90 | 11.1132 | 0.6131 | 0.4892 | 0.7579 | **9.823** |
| 135 | 0.95 | 11.7306 | 0.6473 | 0.5164 | 0.7755 | **11.133** |
| 135 | 1.00 | 12.3480 | 0.6683 | 0.5332 | 0.7740 | **11.785** |
| 135 | 1.05 | 12.9654 | 0.6944 | 0.5541 | 0.7791 | **12.703** |
| 135 | 1.10 | 13.5828 | 0.7233 | 0.5771 | 0.7855 | **13.742** |
| 180 | 0.90 | 11.1132 | 0.5159 | 0.4116 | 0.7772 | **7.152** |
| 180 | 0.95 | 11.7306 | 0.5395 | 0.4304 | 0.7853 | **7.921** |
| 180 | 1.00 | 12.3480 | 0.5677 | 0.4529 | 0.7949 | **8.870** |
| 180 | 1.05 | 12.9654 | 0.5978 | 0.4770 | 0.8045 | **9.908** |
| 180 | 1.10 | 13.5828 | 0.6304 | 0.5030 | 0.8145 | **11.057** |

**KCS_L0p950 (Group 10 variant, Lpp = 218.5 m)**

| β (°) | Fn/Fn_des | U (m/s) | σ_a (m/s²) | ā (m/s²) | ω̄_e (rad/s) | MSI (%) |
|:--|:--|:--|:--|:--|:--|:--|
| 90  | 0.90 | 10.8318 | 0.4318 | 0.3445 | 0.6658 | **3.842** |
| 90  | 0.95 | 11.4336 | 0.4331 | 0.3455 | 0.6658 | **3.869** |
| 90  | 1.00 | 12.0353 | 0.4343 | 0.3465 | 0.6658 | **3.896** |
| 90  | 1.05 | 12.6371 | 0.4355 | 0.3475 | 0.6660 | **3.923** |
| 90  | 1.10 | 13.2389 | 0.4368 | 0.3485 | 0.6661 | **3.951** |
| 135 | 0.90 | 10.8318 | 0.6305 | 0.5031 | 0.7527 | **10.285** |
| 135 | 0.95 | 11.4336 | 0.6512 | 0.5196 | 0.7608 | **11.050** |
| 135 | 1.00 | 12.0353 | 0.6913 | 0.5516 | 0.7691 | **12.452** |
| 135 | 1.05 | 12.6371 | 0.7205 | 0.5749 | 0.7746 | **13.486** |
| 135 | 1.10 | 13.2389 | 0.7462 | 0.5954 | 0.7793 | **14.403** |
| 180 | 0.90 | 10.8318 | 0.5548 | 0.4427 | 0.7724 | **8.241** |
| 180 | 0.95 | 11.4336 | 0.5823 | 0.4646 | 0.7818 | **9.179** |
| 180 | 1.00 | 12.0353 | 0.6090 | 0.4859 | 0.7919 | **10.126** |
| 180 | 1.05 | 12.6371 | 0.6369 | 0.5082 | 0.8005 | **11.115** |
| 180 | 1.10 | 13.2389 | 0.6675 | 0.5326 | 0.8082 | **12.200** |

> **Range across all 30 cases: 3.84 % (L0p950, β = 90°, Fn = 0.90) → 14.40 % (L0p950, β = 135°, Fn = 1.10).**

---

### Deliverable 4 — Figures

The four required plots are listed in the [figures table above](#-figures-in-this-folder-figures). Higher-resolution versions with consistent styling (PDF + PNG, 300 dpi) are in [`Extra Research/figures/`](../Extra%20Research/figures/).

---

### Deliverable 5 — Discussion

*(Report §7; key findings summarised here for reference)*

**5a — Hull comparison**
The shorter KCS_L0p950 is more uncomfortable in head and bow-quartering seas (MSI higher by ≈ 1–1.4 pp at Fn_des) but marginally better in beam seas (≈ 0.1 pp lower). The peak wavelength λ_p = 2πg/ω_p² ≈ 271 m is 1.24 L_pp for the variant vs 1.18 L_pp for the baseline — the shorter hull sits closer to the most pitch-exciting wavelength.

**5b — Speed effect**
At 135°/180° MSI rises monotonically with speed because increasing U shifts ω̄_e upward toward ≈ 1 rad/s (the most provocative frequency in the O'Hanlon–McCauley model, where μ_MSI is minimised). At 90° ω_e = ω (speed-independent), so MSI varies by < 0.2 pp across all speeds.

**5c — Heading effect**
135°/180° cases give 1.8 – 3.7× the beam-sea MSI. At the aft cabin (x_P = −0.35 L_pp), the pitch term x_P ξ₅ dominates the vertical acceleration; beam seas excite mainly heave, which is smaller. Bow-quartering often exceeds head seas here because roll and pitch cross-terms add constructively.

**5d — Cabin position** *(KCS_L0p950, Fn_des)*

| Study | β | x_P (m) | y_P (m) | ζ_v | MSI (%) |
|:--|:--|:--|:--|:--|:--|
| `cabin_midship` | 180° | 0 | +13.2146 | 0.05 | **6.410** |
| `cabin_centreline` | 90° | −76.475 | 0 | 0.05 | **3.727** |
| `cabin_starboard` | 90° | −76.475 | −13.2146 | 0.05 | **3.574** |
| nominal | 90° | −76.475 | +13.2146 | 0.05 | **3.896** |

Moving to midship cuts head-sea MSI from 10.1% to 6.4% — pitch dominates aft cabins. In beam seas, the port cabin (3.90%) is worse than centreline (3.73%) or starboard (3.57%) because waves arrive from starboard and the roll term adds to heave on the port side. **Recommended cabin location: near midship, on the starboard side, as low as structurally permitted.**

**5e — Roll damping** *(KCS_L0p950, β = 90°, Fn_des)*

| Study | ζ_v | MSI (%) |
|:--|:--|:--|
| `damping_low` | 0.025 | **3.893** |
| nominal | 0.050 | **3.896** |
| `damping_high` | 0.100 | **3.901** |

Tripling ζ_v changes MSI by < 0.01 pp. Roll resonance at ω_n ≈ 0.165 rad/s (T_n ≈ 38 s) lies far below the wave spectrum where S(ω_n) ≈ 0, so roll damping has negligible effect on vertical-acceleration MSI.

---

## 🔢 Key Numerical Results at a Glance

| Quantity | KCS_L1p000 | KCS_L0p950 |
|:--|:--|:--|
| L_pp (m) | 230.0 | 218.5 |
| B (m) | 32.20 | 33.04 |
| ω_n (rad/s) | 0.1646 | 0.1626 |
| B_v,44 (N·m·s) | 1.819 × 10⁸ | 1.889 × 10⁸ |
| MSI range, β = 90° (%) | 3.95 – 4.05 | 3.84 – 3.95 |
| MSI range, β = 135° (%) | 9.82 – 13.74 | 10.28 – 14.40 |
| MSI range, β = 180° (%) | 7.15 – 11.06 | 8.24 – 12.20 |

---

## ⚙️ How to Reproduce the Report

```bash
# 1. Install dependencies
python -m pip install -r requirements.txt

# 2. Download HydRA data (if not already in data/raw/)
python tools/get_data.py KCS_L0p950

# 3. Compute all results
python src/msi.py KCS_L0p950

# 4. Build report PDF
python src/build_report.py

# 5. Verify format (should score 20/20)
python tools/check_submission.py
```

Or, to regenerate everything including extra tables and figures:

```bash
python "Extra Research/run_all.py"
```

---

## 📚 Sea State (JONSWAP, IACS Rec. 34 North Atlantic)

| Parameter | Value |
|:--|:--|
| H_s | 4.5 m |
| T₀_m₋₁ (energy period) | 11.5 s |
| T_p (peak period) | 13.184 s |
| ω_p | 0.47659 rad/s |
| γ (JONSWAP peak factor) | 1.5 |
| α (normalisation) | 0.88633 |
| m₀ (variance) | 1.2653 m² |
| 2π m₋₁/m₀ (check) | 11.501 s ≈ 11.5 ✓ |

---

## ⚠️ Placeholders (to fill before submission)

The report title page contains the following blanks — the group must complete them:

- `<<NAME 1>>` through `<<NAME 4>>` — group members' full names
- `<<ROLL NO>>` — corresponding roll numbers
- `<<SUBMISSION DATE>>` — final submission date
- Appendix E: Tools-used statement (AI assistance disclosure per course policy)
