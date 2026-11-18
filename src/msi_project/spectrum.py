"""IACS-normalised JONSWAP spectrum and spectral moments."""
from __future__ import annotations

import numpy as np

from .constants import GAMMA, HS, T0M1


def jonswap_params(Hs: float = HS, T0m1: float = T0M1, gamma: float = GAMMA) -> dict[str, float]:
    """Return Tp, ωp and α for the project JONSWAP normalisation."""
    Tp = T0m1 / (0.7757 + 0.0965 * np.sqrt(gamma) - 0.0144 * gamma)
    omega_p = 2.0 * np.pi / Tp
    alpha = 0.2 / (0.1160 + 0.0594 * np.sqrt(gamma) + 0.0246 * gamma)
    return {"Tp": float(Tp), "omega_p": float(omega_p), "alpha": float(alpha),
            "Hs": float(Hs), "T0m1": float(T0m1), "gamma": float(gamma)}


def jonswap(omega, Hs: float = HS, T0m1: float = T0M1, gamma: float = GAMMA) -> np.ndarray:
    """One-sided JONSWAP S(ω) [m² s] on the wave-frequency axis."""
    o = np.asarray(omega, dtype=np.float64)
    p = jonswap_params(Hs, T0m1, gamma)
    wp = p["omega_p"]
    alpha = p["alpha"]
    safe = np.maximum(o, 1e-30)
    sigma = np.where(o <= wp, 0.07, 0.09)
    peak = np.exp(-0.5 * ((safe / wp - 1.0) / sigma) ** 2)
    S = (
        alpha
        * (5.0 / 16.0)
        * Hs**2
        * wp**4
        * safe ** (-5)
        * np.exp(-1.25 * (wp / safe) ** 4)
        * gamma**peak
    )
    return np.where(o > 0.0, S, 0.0)


def spectral_moments(omega, S, orders=(0, -1)) -> dict[str, float]:
    """Trapezoidal moments m_n = ∫ ω^n S(ω) dω for the requested orders."""
    omega = np.asarray(omega, dtype=np.float64)
    S = np.asarray(S, dtype=np.float64)
    out = {}
    for n in orders:
        if n >= 0:
            integrand = (omega**n) * S
        else:
            integrand = np.zeros_like(S)
            mask = omega > 0
            integrand[mask] = (omega[mask] ** n) * S[mask]
        out[f"m{n}"] = float(np.trapezoid(integrand, omega))
    return out


def check_sea_state(omega_max: float = 6.0, n: int = 60000) -> dict[str, float]:
    """Compute acceptance scalars for Section 5.1."""
    p = jonswap_params()
    omega = np.linspace(1e-6, omega_max, n)
    S = jonswap(omega)
    moms = spectral_moments(omega, S, orders=(0, -1))
    m0 = moms["m0"]
    mm1 = moms["m-1"]
    T0m1_num = 2.0 * np.pi * mm1 / m0
    S_wp = float(jonswap(np.array([p["omega_p"]]))[0])
    return {
        **p,
        "S_wp": S_wp,
        "m0": m0,
        "T0m1_num": T0m1_num,
        "m0_Hs2_16": HS**2 / 16.0,
    }
