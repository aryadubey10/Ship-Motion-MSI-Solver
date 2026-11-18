"""Encounter frequency, 6-DOF frequency-domain solve, roll natural frequency."""
from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator

from .constants import G


def encounter_frequency(omega, U: float, beta_rad: float, g: float = G) -> np.ndarray:
    """Deep-water encounter frequency ωe = ω − (ω² U / g) cos β."""
    o = np.asarray(omega, dtype=np.float64)
    return o - (o**2 * U / g) * np.cos(beta_rad)


def solve_rao(
    omega,
    U: float,
    beta_rad: float,
    M: np.ndarray,
    C: np.ndarray,
    A: np.ndarray,
    B: np.ndarray,
    F: np.ndarray,
    Bv44: float = 0.0,
    g: float = G,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve Z ξ̂ = F̂ on a wave-frequency grid.

    Parameters
    ----------
    A, B : (nω, 6, 6)
    F : (nω, 6) complex
    Returns
    -------
    xi : (nω, 6) complex RAO
    we : (nω,) encounter frequencies
    """
    o = np.asarray(omega, dtype=np.float64)
    we = encounter_frequency(o, U, beta_rad, g=g)
    n = len(o)
    Bd = np.array(B, dtype=np.float64, copy=True)
    Bd[:, 3, 3] += Bv44
    Z = (
        -(we**2)[:, None, None] * (M[None, :, :] + A)
        + 1j * we[:, None, None] * Bd
        + C[None, :, :]
    )
    # Batched solve: Z[i] @ xi[i] = F[i]
    xi = np.linalg.solve(Z, F[..., None])[..., 0]
    assert xi.shape == (n, 6)
    return xi, we


def roll_natural_frequency(
    w_raw: np.ndarray,
    A44_raw: np.ndarray,
    M44: float,
    C44: float,
    tol: float = 1e-10,
    max_iter: int = 100,
) -> tuple[float, float]:
    """Fixed-point iteration for ωn using PCHIP A44(ω).

    PCHIP is used here because ωn lies in 0.1–0.2 rad/s where B44/A44 vary
    steeply; monotone piecewise cubic avoids cubic-spline overshoot.
    The fine-grid hydrodynamic interpolation remains CubicSpline (Section 5.5).

    Returns (ωn, A44(ωn)).
    """
    spl = PchipInterpolator(w_raw, A44_raw, extrapolate=False)
    A = 0.0
    wn = 0.0
    for _ in range(max_iter):
        wn_new = float(np.sqrt(C44 / (M44 + A)))
        if abs(wn_new - wn) < tol:
            wn = wn_new
            break
        wn = wn_new
        if wn < w_raw[0] or wn > w_raw[-1]:
            raise ValueError(f"ωn={wn} outside data range [{w_raw[0]}, {w_raw[-1]}]")
        A = float(spl(wn))
    else:
        raise RuntimeError("roll natural frequency iteration did not converge")
    return wn, A


def roll_damping_from_pack(
    pack: dict,
    zeta_v: float = 0.05,
) -> dict[str, float]:
    """Compute ωn, B_cr, B_v,44 and ζ_wave for one hull.

    Uses A44 and B44 in beam seas (β=90°) at design speed (Fn/Fn_des=1.00).
    """
    ib = pack["beta_index"][90]
    iu = pack["speed_index"][1.0]
    w = pack["w"]
    A44 = pack["AM"][:, 3, 3, ib, iu]
    B44 = pack["BD"][:, 3, 3, ib, iu]
    M44 = float(pack["M"][3, 3])
    C44 = float(pack["C"][3, 3])

    wn, A_wn = roll_natural_frequency(w, A44, M44, C44)
    B_cr = 2.0 * np.sqrt(C44 * (M44 + A_wn))
    Bv44 = zeta_v * B_cr
    # PCHIP for B44(ωn): cubic spline overshoots between 0.1 and 0.2 rad/s.
    B_wave = float(PchipInterpolator(w, B44, extrapolate=False)(wn))
    zeta_wave = B_wave / B_cr
    return {
        "hull": pack["member"],
        "omega_n": float(wn),
        "A44_wn": float(A_wn),
        "B_cr": float(B_cr),
        "Bv44": float(Bv44),
        "zeta_v": float(zeta_v),
        "zeta_wave": float(zeta_wave),
        "B44_wn": float(B_wave),
    }


def reproduce_hydra(pack: dict, rtol: float = 1e-6) -> dict[str, float]:
    """Solve on the raw grid with Bv=0 and compare to HydRA RAOs."""
    max_rel = 0.0
    worst = None
    for beta_deg, ib in pack["beta_index"].items():
        beta_rad = float(pack["beta"][ib])
        for ratio, iu in pack["speed_index"].items():
            U = float(pack["speed"][iu])
            A = pack["AM"][:, :, :, ib, iu]
            B = pack["BD"][:, :, :, ib, iu]
            F = pack["F"][:, :, ib, iu]
            xi, _ = solve_rao(pack["w"], U, beta_rad, pack["M"], pack["C"], A, B, F, Bv44=0.0)
            ref = pack["rao"][:, :, ib, iu]
            for mode in range(6):
                denom = np.max(np.abs(ref[:, mode]))
                if denom < 1e-30:
                    err = float(np.max(np.abs(xi[:, mode] - ref[:, mode])))
                    rel = err  # absolute when HydRA is ~0
                else:
                    rel = float(np.max(np.abs(xi[:, mode] - ref[:, mode])) / denom)
                if rel > max_rel:
                    max_rel = rel
                    worst = (pack["member"], beta_deg, ratio, mode, rel)
    return {"max_rel_error": max_rel, "worst": worst}
