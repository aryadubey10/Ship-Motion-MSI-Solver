"""Deliverable figures (PDF + 300 dpi PNG) and wave-spectrum plot."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm

from .cabin import evaluate_case
from .clean import short_wave_taper, taper_omega
from .constants import FIGURES, FN_RATIOS, G, HEADINGS_DEG
from .motions import solve_rao
from .msi import msi_percent
from .spectrum import jonswap, jonswap_params

# Baseline solid, variant dashed; same colour per heading.
HEADING_COLORS = {90: "#1f4e79", 135: "#b85c38", 180: "#2f6b3c"}
HULL_LS = {"KCS_L1p000": "-", "KCS_L0p950": "--"}


def _save(fig: plt.Figure, stem: str, out_dir: Path = FIGURES) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig_wave_spectrum(out_dir: Path = FIGURES) -> Path:
    p = jonswap_params()
    w = np.linspace(0.01, 2.5, 2500)
    S = jonswap(w)
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(w, S, color="#1f4e79", lw=1.8, label="JONSWAP γ=1.5")
    ax.axvline(p["omega_p"], color="0.4", ls=":", lw=1, label=f"ωp={p['omega_p']:.4f} rad/s")
    ax.set_xlim(0, 2.0)
    ax.set_xlabel("Wave frequency ω (rad/s)")
    ax.set_ylabel("S(ω) (m²·s)")
    ax.set_title("Design sea state: Hs = 4.5 m, T₀₋₁ = 11.5 s")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    _save(fig, "fig_wave_spectrum", out_dir)
    return out_dir / "fig_wave_spectrum.pdf"


def fig1_roll_rao_beam(packs: dict, roll_df: pd.DataFrame, out_dir: Path = FIGURES) -> Path:
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for hull, pack in packs.items():
        ib = pack["beta_index"][90]
        iu = pack["speed_index"][1.0]
        beta_rad = float(pack["beta"][ib])
        U = float(pack["speed"][iu])
        A = pack["AM"][:, :, :, ib, iu]
        B = pack["BD"][:, :, :, ib, iu]
        F = pack["F"][:, :, ib, iu]
        Bv = float(roll_df.loc[roll_df["hull"] == hull, "Bv44"].iloc[0])
        wn = float(roll_df.loc[roll_df["hull"] == hull, "omega_n"].iloc[0])
        for Bv44, label_extra, alpha in ((0.0, "no Bv", 0.55), (Bv, "with Bv", 1.0)):
            xi, we = solve_rao(pack["w"], U, beta_rad, pack["M"], pack["C"], A, B, F, Bv44=Bv44)
            ax.semilogy(
                we,
                np.abs(xi[:, 3]),
                ls=HULL_LS[hull],
                color="#1f4e79" if Bv44 == 0 else "#b85c38",
                alpha=alpha,
                lw=1.5,
                label=f"{hull}, {label_extra}",
            )
        ax.axvline(wn, color="0.5", ls=":", lw=0.9)
        ax.text(wn, ax.get_ylim()[0] * 1.5 if False else 1e-3, "", fontsize=8)
    # Mark ωn for baseline
    wn0 = float(roll_df.loc[roll_df["hull"] == "KCS_L1p000", "omega_n"].iloc[0])
    ax.axvline(wn0, color="0.35", ls=":", lw=1.0, label=f"ωn (L1p000)={wn0:.3f}")
    ax.set_xlabel("Encounter frequency ωe (rad/s)")
    ax.set_ylabel("|ξ₄| (rad/m)")
    ax.set_title("Roll RAO, beam seas, Fn/Fn_des = 1.00")
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    _save(fig, "fig1_roll_rao_beam", out_dir)
    return out_dir / "fig1_roll_rao_beam.pdf"


def fig2_msi_vs_speed(main_df: pd.DataFrame, out_dir: Path = FIGURES) -> Path:
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for hull in ("KCS_L1p000", "KCS_L0p950"):
        for beta in HEADINGS_DEG:
            sub = main_df[(main_df["hull"] == hull) & (main_df["beta_deg"] == beta)].sort_values(
                "fn_ratio"
            )
            ax.plot(
                sub["fn_ratio"],
                sub["msi_pct"],
                ls=HULL_LS[hull],
                color=HEADING_COLORS[beta],
                marker="o",
                ms=4,
                lw=1.6,
                label=f"{hull}, {beta}°",
            )
    ax.set_xlabel("Fn / Fn_des")
    ax.set_ylabel("MSI (%)")
    ax.set_title("MSI versus speed")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    _save(fig, "fig2_msi_vs_speed", out_dir)
    return out_dir / "fig2_msi_vs_speed.pdf"


def fig3_accel_spectra(packs: dict, roll_df: pd.DataFrame, out_dir: Path = FIGURES) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), sharey=False)
    for ax, beta in zip(axes, HEADINGS_DEG):
        for hull, pack in packs.items():
            Bv = float(roll_df.loc[roll_df["hull"] == hull, "Bv44"].iloc[0])
            out = evaluate_case(pack, beta, 1.0, Bv, apply_taper=True)
            ax.plot(
                out["we_plot"],
                out["S_aa"],
                ls=HULL_LS[hull],
                color=HEADING_COLORS[beta],
                lw=1.5,
                label=hull,
            )
        ax.set_xlabel("ωe (rad/s)")
        ax.set_ylabel(r"$S_{aa}(\omega_e)$ [(m/s²)²·s]")
        ax.set_title(f"β = {beta}°")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7)
    fig.suptitle("Cabin acceleration spectra at Fn/Fn_des = 1.00", y=1.02)
    fig.tight_layout()
    _save(fig, "fig3_accel_spectra", out_dir)
    return out_dir / "fig3_accel_spectra.pdf"


def fig4_msi_plane(main_df: pd.DataFrame, out_dir: Path = FIGURES) -> Path:
    we = np.linspace(0.5, 1.2, 200)
    aa = np.linspace(0.2, 0.8, 200)
    WE, AA = np.meshgrid(we, aa)
    # Vectorised MSI
    mu = -0.819 + 2.32 * (np.log10(WE)) ** 2
    Z = 100.0 * norm.cdf((np.log10(AA / G) - mu) / 0.4)

    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    cs = ax.contour(WE, AA, Z, levels=[2, 5, 10, 15, 20, 30], colors="0.55", linewidths=0.9)
    ax.clabel(cs, inline=True, fontsize=8, fmt="%g%%")

    markers = {90: "o", 135: "s", 180: "^"}
    fills = {"KCS_L1p000": "full", "KCS_L0p950": "none"}
    for hull in ("KCS_L1p000", "KCS_L0p950"):
        for beta in HEADINGS_DEG:
            sub = main_df[(main_df["hull"] == hull) & (main_df["beta_deg"] == beta)]
            sizes = 40 + 80 * (sub["fn_ratio"] - 0.9) / 0.2
            ax.scatter(
                sub["we_bar"],
                sub["abar"],
                s=sizes,
                marker=markers[beta],
                facecolors=HEADING_COLORS[beta] if hull == "KCS_L1p000" else "none",
                edgecolors=HEADING_COLORS[beta],
                linewidths=1.3,
                label=f"{hull}, {beta}°",
            )
    ax.set_xlim(0.5, 1.2)
    ax.set_ylim(0.2, 0.8)
    ax.set_xlabel(r"Mean encounter frequency $\bar{\omega}_e$ (rad/s)")
    ax.set_ylabel(r"Mean absolute acceleration $\bar{a}$ (m/s²)")
    ax.set_title("30 cases on the MSI plane (marker size ∝ speed)")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=7, ncol=2, loc="upper left")
    fig.tight_layout()
    _save(fig, "fig4_msi_plane", out_dir)
    return out_dir / "fig4_msi_plane.pdf"


def make_all_figures(packs: dict, roll_df: pd.DataFrame, main_df: pd.DataFrame) -> list[Path]:
    paths = [
        fig_wave_spectrum(),
        fig1_roll_rao_beam(packs, roll_df),
        fig2_msi_vs_speed(main_df),
        fig3_accel_spectra(packs, roll_df),
        fig4_msi_plane(main_df),
    ]
    return paths
