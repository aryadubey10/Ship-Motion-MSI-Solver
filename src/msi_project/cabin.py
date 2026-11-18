"""Cabin vertical acceleration RAO, response spectrum and moments."""
from __future__ import annotations

import numpy as np

from .clean import short_wave_taper
from .constants import G
from .msi import from_sigma
from .motions import encounter_frequency, solve_rao
from .spectrum import jonswap


def cabin_vertical_displacement(xi: np.ndarray, x_P: float, y_P: float) -> np.ndarray:
    """ŝ3 = ξ̂3 + y_P ξ̂4 − x_P ξ̂5 (per metre wave amplitude)."""
    return xi[:, 2] + y_P * xi[:, 3] - x_P * xi[:, 4]


def cabin_acceleration_rao(
    xi: np.ndarray,
    we: np.ndarray,
    x_P: float,
    y_P: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Ĥa = −ωe² ŝ3. Returns (Ha, s3)."""
    s3 = cabin_vertical_displacement(xi, x_P, y_P)
    Ha = -(we**2) * s3
    return Ha, s3


def acceleration_moments(
    omega: np.ndarray,
    Ha: np.ndarray,
    we: np.ndarray,
    S_wave: np.ndarray | None = None,
) -> dict[str, float]:
    """m_n^a = ∫ |ωe|^n |Ĥa|² S(ω) dω over the supplied ω grid (trapezoid)."""
    if S_wave is None:
        S_wave = jonswap(omega)
    Sa = (np.abs(Ha) ** 2) * S_wave
    m0 = float(np.trapezoid(Sa, omega))
    m1 = float(np.trapezoid(np.abs(we) * Sa, omega))
    sigma_a = float(np.sqrt(max(m0, 0.0)))
    we_bar = float(m1 / m0) if m0 > 0 else float("nan")
    return {"m0": m0, "m1": m1, "sigma_a": sigma_a, "we_bar": we_bar, "S_a_wave": Sa}


def encounter_spectrum(
    omega: np.ndarray,
    Ha: np.ndarray,
    U: float,
    beta_rad: float,
    S_wave: np.ndarray | None = None,
    g: float = G,
) -> tuple[np.ndarray, np.ndarray, float]:
    """S_aa(ωe) = |Ĥa|² S(ω) / |dωe/dω|; returns (ωe, S_aa, ∫S_aa dωe)."""
    if S_wave is None:
        S_wave = jonswap(omega)
    we = encounter_frequency(omega, U, beta_rad, g=g)
    dwe = 1.0 - 2.0 * omega * U * np.cos(beta_rad) / g
    Sa_wave = (np.abs(Ha) ** 2) * S_wave
    Saa = Sa_wave / np.maximum(np.abs(dwe), 1e-30)
    # Integrate vs encounter frequency (monotonic for β∈{90,135,180}).
    order = np.argsort(we)
    m0_e = float(np.trapezoid(Saa[order], we[order]))
    return we, Saa, m0_e


def evaluate_case(
    pack: dict,
    beta_deg: int,
    fn_ratio: float,
    Bv44: float,
    x_P: float | None = None,
    y_P: float | None = None,
    apply_taper: bool = True,
) -> dict:
    """Full cabin acceleration + MSI evaluation for one (hull, β, U) case."""
    ib = pack["beta_index"][beta_deg]
    iu = pack["speed_index"][float(fn_ratio)]
    beta_rad = float(pack["beta"][ib])
    U = float(pack["speed"][iu])
    x_P = pack["x_P"] if x_P is None else x_P
    y_P = pack["y_P"] if y_P is None else y_P

    A = pack["AM"][:, :, :, ib, iu]
    B = pack["BD"][:, :, :, ib, iu]
    F = pack["F"][:, :, ib, iu]
    xi, we = solve_rao(pack["w"], U, beta_rad, pack["M"], pack["C"], A, B, F, Bv44=Bv44)
    Ha, s3 = cabin_acceleration_rao(xi, we, x_P, y_P)

    if apply_taper and beta_deg in (135, 180):
        omega_t = pack.get("omega_t")
        if omega_t is None:
            from .clean import taper_omega

            omega_t = taper_omega(pack["Lpp"])
        Ha = Ha * short_wave_taper(pack["w"], omega_t)

    moms = acceleration_moments(pack["w"], Ha, we)
    msi = from_sigma(moms["sigma_a"], moms["we_bar"])
    we_plot, Saa, m0_e = encounter_spectrum(pack["w"], Ha, U, beta_rad)
    return {
        "hull": pack["member"],
        "beta_deg": beta_deg,
        "fn_ratio": float(fn_ratio),
        "U_mps": U,
        "x_P_m": float(x_P),
        "y_P_m": float(y_P),
        "Bv44": float(Bv44),
        "xi": xi,
        "we": we,
        "Ha": Ha,
        "s3": s3,
        "S_aa": Saa,
        "we_plot": we_plot,
        "m0_encounter": m0_e,
        **moms,
        **msi,
    }


def contribution_split(
    pack: dict,
    beta_deg: int,
    fn_ratio: float,
    Bv44: float,
    x_P: float | None = None,
    y_P: float | None = None,
    apply_taper: bool = True,
) -> dict:
    """m0 of heave, roll, pitch terms of Ĥa alone, plus pairwise cross-terms."""
    ib = pack["beta_index"][beta_deg]
    iu = pack["speed_index"][float(fn_ratio)]
    beta_rad = float(pack["beta"][ib])
    U = float(pack["speed"][iu])
    x_P = pack["x_P"] if x_P is None else x_P
    y_P = pack["y_P"] if y_P is None else y_P

    A = pack["AM"][:, :, :, ib, iu]
    B = pack["BD"][:, :, :, ib, iu]
    F = pack["F"][:, :, ib, iu]
    xi, we = solve_rao(pack["w"], U, beta_rad, pack["M"], pack["C"], A, B, F, Bv44=Bv44)

    # Acceleration contributions: Ha = −ωe² (ξ3 + y ξ4 − x ξ5)
    h = -(we**2) * xi[:, 2]
    r = -(we**2) * (y_P * xi[:, 3])
    p = -(we**2) * (-x_P * xi[:, 4])
    if apply_taper and beta_deg in (135, 180):
        from .clean import taper_omega

        omega_t = pack.get("omega_t", taper_omega(pack["Lpp"]))
        tap = short_wave_taper(pack["w"], omega_t)
        h, r, p = h * tap, r * tap, p * tap

    S = jonswap(pack["w"])

    def m0_of(a, b=None):
        if b is None:
            return float(np.trapezoid(np.abs(a) ** 2 * S, pack["w"]))
        return float(np.trapezoid(np.real(a * np.conj(b)) * S, pack["w"]))

    return {
        "hull": pack["member"],
        "beta_deg": beta_deg,
        "fn_ratio": float(fn_ratio),
        "m0_heave": m0_of(h),
        "m0_roll": m0_of(r),
        "m0_pitch": m0_of(p),
        "m0_heave_roll": 2 * m0_of(h, r),
        "m0_heave_pitch": 2 * m0_of(h, p),
        "m0_roll_pitch": 2 * m0_of(r, p),
        "m0_total": m0_of(h + r + p),
    }
