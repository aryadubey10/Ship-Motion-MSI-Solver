"""Case loops: 30 main cases, sensitivities, taper effect, contributions."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .cabin import contribution_split, evaluate_case
from .clean import refine_pack, write_diagnostic_overlays, write_spikes_csv
from .constants import (
    COLS_MSI_CASES,
    COLS_ROLL,
    COLS_SENSITIVITY,
    EXTRA_RESULTS,
    FN_RATIOS,
    HEADINGS_DEG,
    MEMBERS,
    RESULTS,
    ZETA_V_BASE,
)
from .io_hydra import load_npz, process_all
from .motions import roll_damping_from_pack


def prepare_packs(force_reload: bool = False) -> dict[str, dict]:
    """Load processed NPZ (building from raw if needed) and refine to fine grid."""
    need = force_reload
    for m in MEMBERS:
        from .constants import DATA_PROCESSED

        if not (DATA_PROCESSED / f"{m}.npz").exists():
            need = True
    if need:
        process_all()
    packs = {}
    for m in MEMBERS:
        raw = load_npz(m)
        # Attach beta/speed index from load_npz
        fine = refine_pack(raw)
        write_spikes_csv(m, fine["spike_records"])
        write_diagnostic_overlays(fine)
        packs[m] = fine
    _write_official_spikes(packs)
    return packs


def _write_official_spikes(packs: dict[str, dict]) -> None:
    """results/spikes.csv in the course template format (one row per replaced ω)."""
    seen: set[tuple] = set()
    rows = []
    for pack in packs.values():
        speeds = list(pack["speed"])
        for rec in pack.get("spike_records", []):
            U = float(rec["U"])
            iu = min(range(len(speeds)), key=lambda i: abs(speeds[i] - U))
            fn = float(FN_RATIOS[iu])
            key = (rec["member"], int(rec["beta_deg"]), fn, round(float(rec["omega"]), 6))
            if key in seen:
                continue
            seen.add(key)
            rows.append(
                {
                    "hull": rec["member"],
                    "beta_deg": int(rec["beta_deg"]),
                    "fn_ratio": fn,
                    "omega": float(rec["omega"]),
                }
            )
    RESULTS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["hull", "beta_deg", "fn_ratio", "omega"]).to_csv(
        RESULTS / "spikes.csv", index=False
    )


def compute_roll_damping(packs: dict[str, dict], zeta_v: float = ZETA_V_BASE) -> pd.DataFrame:
    rows = []
    for m, pack in packs.items():
        # Roll damping uses raw-grid A44 (before fine spline is fine — refine keeps w_raw).
        raw_like = {
            **pack,
            "w": pack["w_raw"],
            "AM": pack["AM_raw"],
            "BD": pack["BD_raw"],
        }
        row = roll_damping_from_pack(raw_like, zeta_v=zeta_v)
        rows.append(row)
    df = pd.DataFrame(rows)
    RESULTS.mkdir(parents=True, exist_ok=True)
    EXTRA_RESULTS.mkdir(parents=True, exist_ok=True)
    df.to_csv(EXTRA_RESULTS / "roll_damping_full.csv", index=False)
    df[COLS_ROLL].to_csv(RESULTS / "roll_damping.csv", index=False)
    return df


def run_main_cases(
    packs: dict[str, dict],
    roll_df: pd.DataFrame,
    zeta_v: float = ZETA_V_BASE,
) -> pd.DataFrame:
    bv = {r["hull"]: float(r["Bv44"]) for _, r in roll_df.iterrows()}
    # Bv44 at baseline ζv; scale if needed
    rows = []
    for m, pack in packs.items():
        Bv44 = bv[m] * (zeta_v / ZETA_V_BASE)
        for beta in HEADINGS_DEG:
            for ratio in FN_RATIOS:
                out = evaluate_case(
                    pack,
                    beta_deg=beta,
                    fn_ratio=ratio,
                    Bv44=Bv44,
                    apply_taper=True,
                )
                rows.append(
                    {
                        "hull": out["hull"],
                        "beta_deg": int(out["beta_deg"]),
                        "fn_ratio": float(ratio),
                        "U_mps": out["U_mps"],
                        "x_P_m": out["x_P_m"],
                        "y_P_m": out["y_P_m"],
                        "zeta_v": zeta_v,
                        "sigma_a": out["sigma_a"],
                        "abar": out["abar"],
                        "we_bar": out["we_bar"],
                        "mu_msi": out["mu_msi"],
                        "msi_pct": out["msi_pct"],
                    }
                )
    df = pd.DataFrame(rows)
    EXTRA_RESULTS.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    df.to_csv(EXTRA_RESULTS / "main_results.csv", index=False)
    df[COLS_MSI_CASES].to_csv(RESULTS / "msi_cases.csv", index=False)
    return df


def run_sensitivity(packs: dict[str, dict], roll_df: pd.DataFrame) -> pd.DataFrame:
    """Sensitivity on KCS_L0p950 at Fn/Fn_des = 1.00."""
    pack = packs["KCS_L0p950"]
    B_cr_row = roll_df.loc[roll_df["hull"] == "KCS_L0p950"].iloc[0]
    B_cr = float(B_cr_row["Bv44"]) / ZETA_V_BASE  # Bv44 = ζv * B_cr
    B = pack["B"]
    Lpp = pack["Lpp"]
    x_aft = -0.35 * Lpp
    y_port = 0.4 * B

    specs = [
        ("cabin_midship", 180, 0.0, y_port, 0.05),
        ("cabin_centreline", 90, x_aft, 0.0, 0.05),
        ("cabin_starboard", 90, x_aft, -y_port, 0.05),
        ("damping_low", 90, x_aft, y_port, 0.025),
        ("damping_high", 90, x_aft, y_port, 0.10),
    ]
    rows = []
    for item, beta, x_P, y_P, zeta in specs:
        Bv44 = zeta * B_cr
        out = evaluate_case(
            pack,
            beta_deg=beta,
            fn_ratio=1.0,
            Bv44=Bv44,
            x_P=x_P,
            y_P=y_P,
            apply_taper=True,
        )
        rows.append(
            {
                "study": item,
                "hull": out["hull"],
                "beta_deg": int(out["beta_deg"]),
                "fn_ratio": 1.00,
                "U_mps": out["U_mps"],
                "x_P_m": out["x_P_m"],
                "y_P_m": out["y_P_m"],
                "zeta_v": zeta,
                "sigma_a": out["sigma_a"],
                "abar": out["abar"],
                "we_bar": out["we_bar"],
                "mu_msi": out["mu_msi"],
                "msi_pct": out["msi_pct"],
            }
        )
    df = pd.DataFrame(rows)
    EXTRA_RESULTS.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    df.to_csv(EXTRA_RESULTS / "sensitivity_full.csv", index=False)
    df[COLS_SENSITIVITY].to_csv(RESULTS / "sensitivity.csv", index=False)
    return df


def run_taper_effect(packs: dict[str, dict], roll_df: pd.DataFrame) -> pd.DataFrame:
    bv = {r["hull"]: float(r["Bv44"]) for _, r in roll_df.iterrows()}
    rows = []
    for m, pack in packs.items():
        Bv44 = bv[m]
        for beta in (135, 180):
            for ratio in FN_RATIOS:
                on = evaluate_case(pack, beta, ratio, Bv44, apply_taper=True)
                off = evaluate_case(pack, beta, ratio, Bv44, apply_taper=False)
                rows.append(
                    {
                        "hull": m,
                        "beta_deg": beta,
                        "fn_ratio": float(ratio),
                        "msi_taper_on": on["msi_pct"],
                        "msi_taper_off": off["msi_pct"],
                        "delta_msi": on["msi_pct"] - off["msi_pct"],
                    }
                )
    df = pd.DataFrame(rows)
    EXTRA_RESULTS.mkdir(parents=True, exist_ok=True)
    df.to_csv(EXTRA_RESULTS / "taper_effect.csv", index=False)
    return df


def run_contributions(packs: dict[str, dict], roll_df: pd.DataFrame) -> pd.DataFrame:
    bv = {r["hull"]: float(r["Bv44"]) for _, r in roll_df.iterrows()}
    rows = []
    for m, pack in packs.items():
        Bv44 = bv[m]
        for beta in HEADINGS_DEG:
            for ratio in FN_RATIOS:
                rows.append(
                    contribution_split(pack, beta, ratio, Bv44, apply_taper=True)
                )
    df = pd.DataFrame(rows)
    EXTRA_RESULTS.mkdir(parents=True, exist_ok=True)
    df.to_csv(EXTRA_RESULTS / "contributions.csv", index=False)
    return df
