"""Fixed physical constants, sea state, hull table and case grid."""
from __future__ import annotations

from pathlib import Path

import numpy as np

# Gravity used by HydRA and the reference MSI values (m/s^2).
G = 9.80665

# Design sea state (IACS Rec. 34 North Atlantic).
HS = 4.5          # significant wave height (m)
T0M1 = 11.5       # energy period (s)
GAMMA = 1.5       # JONSWAP peak enhancement

# Viscous roll damping ratios.
ZETA_V_BASE = 0.05
ZETA_V_SENS = (0.025, 0.10)

# Headings (deg) and Froude-number ratios.
HEADINGS_DEG = (90, 135, 180)
FN_RATIOS = (0.90, 0.95, 1.00, 1.05, 1.10)
FN_DES = 0.260

# Fine wave-frequency grid (rad/s): 0.1 : 0.001 : 1.5 (1401 points).
OMEGA_FINE = np.linspace(0.1, 1.5, 1401)

# Displacement volume from the HydRA family table (m^3); both members share it.
NABLA = 51953.142819580426

# Nominal hull particulars (prefer ZIP / JSON values when they disagree).
HULLS = {
    "KCS_L1p000": {
        "scale": 1.00,
        "Lpp": 230.0,
        "B": 32.20,
        "T": 10.8002323,
        "nabla": NABLA,
        "CB": 0.651,
        "GM_nominal": 0.60,
        "kxx_over_B": 0.40,
        "kyy_over_Lpp": 0.25,
        "Fn_des": FN_DES,
        "x_P": -0.35 * 230.0,   # -80.50 m
        "y_P": 0.4 * 32.20,     # +12.88 m (port)
    },
    "KCS_L0p950": {
        "scale": 0.95,
        "Lpp": 218.5,
        "B": 32.20 / np.sqrt(0.95),   # 33.0365...
        "T": 10.80 / np.sqrt(0.95),   # 11.0806...
        "nabla": NABLA,
        "CB": 0.651,
        "GM_nominal": 0.62,
        "kxx_over_B": 0.40,
        "kyy_over_Lpp": 0.25,
        "Fn_des": FN_DES,
        "x_P": -0.35 * 218.5,         # -76.475 m
        "y_P": 0.4 * (32.20 / np.sqrt(0.95)),  # +13.2146 m (port)
    },
}

# Repository paths. Official autograder files stay in results/; everything else
# that is not in the course template lives under "Extra Research/".
ROOT = Path(__file__).resolve().parents[2]
DATA_RAW = ROOT / "data" / "raw"
DATA_PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
EXTRA_INFO = ROOT / "Extra Research"
EXTRA_RESULTS = EXTRA_INFO / "results"
FIGURES = EXTRA_INFO / "figures"
DIAGNOSTICS = FIGURES / "diagnostics"
REFERENCE_CSV = EXTRA_INFO / "reference_results_group10.csv"

# Template / autograder CSV columns (tools/msi_checks.py).
COLS_MSI_CASES = [
    "hull", "beta_deg", "fn_ratio", "U_mps", "sigma_a", "abar", "we_bar", "mu_msi", "msi_pct",
]
COLS_ROLL = ["hull", "omega_n", "zeta_wave", "Bv44"]
COLS_SENSITIVITY = [
    "study", "hull", "beta_deg", "fn_ratio", "x_P_m", "y_P_m", "zeta_v",
    "sigma_a", "abar", "we_bar", "mu_msi", "msi_pct",
]
COLS_SPIKES = ["hull", "beta_deg", "fn_ratio", "omega"]

MEMBERS = ("KCS_L1p000", "KCS_L0p950")
MODE_NAMES = ("surge", "sway", "heave", "roll", "pitch", "yaw")
