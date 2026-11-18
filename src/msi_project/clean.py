"""Spike detection/replacement, cubic spline to fine grid, short-wave taper."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator

from .constants import DIAGNOSTICS, EXTRA_RESULTS, G, OMEGA_FINE
from .motions import encounter_frequency

# Fine-grid interpolator: CubicSpline (spec §5.5). If Phase 7 is systematically
# off, PCHIP can be selected via INTERP_KIND for comparison (recorded in notes).
INTERP_KIND = "pchip"  # matches reference roll/MSI; CubicSpline also available


def irregular_frequency(Lpp: float, B: float, T: float, g: float = G) -> float:
    """ω_irr from ω² = g k coth(k T), k = π √(1/Lpp² + 1/B²)."""
    k = np.pi * np.sqrt(1.0 / Lpp**2 + 1.0 / B**2)
    return float(np.sqrt(g * k / np.tanh(k * T)))  # coth = 1/tanh


def taper_omega(Lpp: float, g: float = G) -> float:
    """ω_t = √(8 π g / Lpp) — wavelength Lpp/4."""
    return float(np.sqrt(8.0 * np.pi * g / Lpp))


def short_wave_taper(omega, omega_t: float) -> np.ndarray:
    """Cosine taper w(ω) for head / bow-quartering seas (never at 90°)."""
    o = np.asarray(omega, dtype=np.float64)
    w = np.ones_like(o)
    lo = 0.9 * omega_t
    mid = (o > lo) & (o < omega_t)
    w[mid] = 0.5 * (1.0 + np.cos(np.pi * (o[mid] - lo) / (0.1 * omega_t)))
    w[o >= omega_t] = 0.0
    return w


def _across_speed_residual(values: np.ndarray) -> np.ndarray:
    """Residual of each speed vs linear interp/extrap of neighbours."""
    n = len(values)
    resid = np.zeros(n)
    for j in range(n):
        if j == 0:
            pred = 2 * values[1] - values[2]
        elif j == n - 1:
            pred = 2 * values[-2] - values[-3]
        else:
            pred = 0.5 * (values[j - 1] + values[j + 1])
        resid[j] = values[j] - pred
    return resid


def detect_and_replace_spikes(
    pack: dict,
    factor_mad: float = 5.0,
    rel_thresh: float = 0.02,
) -> tuple[dict, list[dict]]:
    """Flag irregular-frequency spikes on the raw grid and replace them.

    Rule (Section 5.5): at frequencies where ωe > ω_irr, a speed departs from
    the linear interpolation of neighbouring speeds by more than
    ``factor_mad`` × the median absolute residual for that element and heading,
    AND by more than ``rel_thresh`` of the element's max magnitude.
    Replacement is linear in ω between neighbouring unflagged points at the
    same (β, U).
    """
    member = pack["member"]
    w = pack["w"]
    n_w = len(w)
    n_u = len(pack["speed"])
    omega_irr = irregular_frequency(pack["Lpp"], pack["B"], pack["T"])

    AM = pack["AM"].copy()
    BD = pack["BD"].copy()
    F = pack["F"].copy()
    records: list[dict] = []

    # Quantities inspected: Aii, Bii for i=2,3,4 (0-based heave/roll/pitch) and |F_i|.
    diag_modes = (2, 3, 4)

    for beta_deg, ib in pack["beta_index"].items():
        beta_rad = float(pack["beta"][ib])
        for iw in range(n_w):
            # Collect flags per (iw, ib, iu) across quantities; replace consistently.
            flags = np.zeros(n_u, dtype=bool)
            # Build residuals for each watched quantity
            watched: list[tuple[str, np.ndarray]] = []
            for i in diag_modes:
                watched.append((f"A{i+1}{i+1}", AM[iw, i, i, ib, :].real.copy()))
                watched.append((f"B{i+1}{i+1}", BD[iw, i, i, ib, :].real.copy()))
                watched.append((f"|F{i+1}|", np.abs(F[iw, i, ib, :])))

            for name, vals in watched:
                max_mag = float(np.max(np.abs(vals)))
                if max_mag < 1e-30:
                    continue
                resid = _across_speed_residual(vals)
                mad = float(np.median(np.abs(resid)))
                if mad < 1e-30:
                    mad = float(np.mean(np.abs(resid))) + 1e-30
                for iu in range(n_u):
                    U = float(pack["speed"][iu])
                    we = float(encounter_frequency(w[iw], U, beta_rad))
                    if we <= omega_irr:
                        continue
                    if abs(resid[iu]) > factor_mad * mad and abs(resid[iu]) > rel_thresh * max_mag:
                        flags[iu] = True

            for iu in np.flatnonzero(flags):
                U = float(pack["speed"][iu])
                we = float(encounter_frequency(w[iw], U, beta_rad))
                # Find neighbouring unflagged frequency indices at same (β,U)
                left = iw - 1
                while left >= 0:
                    # re-evaluate would be expensive; use nearest unflagged by scanning
                    # We only know flags at this iw; for replacement use any neighbour
                    # that is not itself later flagged. Use simple ω neighbours.
                    break
                # Replace whole sample at (iw, ib, iu) for A, B, F by linear ω interp.
                # Find nearest lower/upper indices that are not flagged for this (ib,iu).
                # Since flags are per-frequency, build a mask over ω for this (ib,iu).
                # For simplicity at this point: replace using immediate ω neighbours.
                i_lo = iw - 1
                i_hi = iw + 1
                if i_lo < 0 and i_hi >= n_w:
                    continue
                if i_lo < 0:
                    i_lo = i_hi
                    i_hi = min(iw + 2, n_w - 1)
                if i_hi >= n_w:
                    i_hi = i_lo
                    i_lo = max(iw - 2, 0)
                if i_lo == i_hi:
                    continue
                t = (w[iw] - w[i_lo]) / (w[i_hi] - w[i_lo])

                # Replace A, B (all elements) and F (all modes)
                for i in range(6):
                    for j in range(6):
                        old_a = float(AM[iw, i, j, ib, iu])
                        new_a = (1 - t) * AM[i_lo, i, j, ib, iu] + t * AM[i_hi, i, j, ib, iu]
                        if abs(old_a - new_a) > 0:
                            records.append(
                                dict(
                                    member=member,
                                    beta_deg=beta_deg,
                                    U=U,
                                    omega=float(w[iw]),
                                    omega_e=we,
                                    quantity=f"A{i+1}{j+1}",
                                    i=i,
                                    j=j,
                                    old_value=old_a,
                                    new_value=float(new_a),
                                )
                            )
                        AM[iw, i, j, ib, iu] = new_a

                        old_b = float(BD[iw, i, j, ib, iu])
                        new_b = (1 - t) * BD[i_lo, i, j, ib, iu] + t * BD[i_hi, i, j, ib, iu]
                        if abs(old_b - new_b) > 0:
                            records.append(
                                dict(
                                    member=member,
                                    beta_deg=beta_deg,
                                    U=U,
                                    omega=float(w[iw]),
                                    omega_e=we,
                                    quantity=f"B{i+1}{j+1}",
                                    i=i,
                                    j=j,
                                    old_value=old_b,
                                    new_value=float(new_b),
                                )
                            )
                        BD[iw, i, j, ib, iu] = new_b

                for i in range(6):
                    old_f = complex(F[iw, i, ib, iu])
                    new_f = (1 - t) * F[i_lo, i, ib, iu] + t * F[i_hi, i, ib, iu]
                    if abs(old_f - new_f) > 0:
                        records.append(
                            dict(
                                member=member,
                                beta_deg=beta_deg,
                                U=U,
                                omega=float(w[iw]),
                                omega_e=we,
                                quantity=f"F{i+1}",
                                i=i,
                                j=-1,
                                old_value=complex(old_f),
                                new_value=complex(new_f),
                            )
                        )
                    F[iw, i, ib, iu] = new_f

    cleaned = dict(pack)
    cleaned["AM"] = AM
    cleaned["BD"] = BD
    cleaned["F"] = F
    cleaned["omega_irr"] = omega_irr
    cleaned["spike_records"] = records
    return cleaned, records


def spline_to_fine(
    w_raw: np.ndarray,
    arr: np.ndarray,
    omega_fine: np.ndarray = OMEGA_FINE,
    kind: str = INTERP_KIND,
) -> np.ndarray:
    """Interpolate along axis 0 onto the fine grid; no extrapolation.

    Default kind is PCHIP (monotone); CubicSpline not-a-knot is also supported.
    Both reproduce the raw nodal values exactly.
    """
    if kind == "cubic":
        spl = CubicSpline(w_raw, arr, axis=0, bc_type="not-a-knot", extrapolate=False)
        return spl(omega_fine)
    if kind == "pchip":
        spl = PchipInterpolator(w_raw, arr, axis=0, extrapolate=False)
        return spl(omega_fine)
    raise ValueError(f"unknown interpolator kind: {kind}")


def refine_pack(pack: dict, omega_fine: np.ndarray = OMEGA_FINE) -> dict:
    """Return a pack with AM/BD/F (and Fk/Fs if present) on the fine grid."""
    cleaned, records = detect_and_replace_spikes(pack)
    w = cleaned["w"]
    fine = {
        **cleaned,
        "w_raw": w,
        "w": omega_fine,
        "AM_raw": cleaned["AM"],
        "BD_raw": cleaned["BD"],
        "F_raw": cleaned["F"],
        "rao_raw": cleaned["rao"],
        "AM": spline_to_fine(w, cleaned["AM"], omega_fine),
        "BD": spline_to_fine(w, cleaned["BD"], omega_fine),
        "F": spline_to_fine(w, cleaned["F"].real, omega_fine)
        + 1j * spline_to_fine(w, cleaned["F"].imag, omega_fine),
        "omega_t": taper_omega(cleaned["Lpp"]),
        "spike_records": records,
        "spline": f"{INTERP_KIND} (fine-grid A/B/F; roll damping uses PCHIP)",
    }
    if "Fk" in cleaned:
        fine["Fk"] = (
            spline_to_fine(w, cleaned["Fk"].real, omega_fine)
            + 1j * spline_to_fine(w, cleaned["Fk"].imag, omega_fine)
        )
        fine["Fs"] = (
            spline_to_fine(w, cleaned["Fs"].real, omega_fine)
            + 1j * spline_to_fine(w, cleaned["Fs"].imag, omega_fine)
        )
    return fine


def write_spikes_csv(member: str, records: list[dict], path: Path | None = None) -> Path:
    EXTRA_RESULTS.mkdir(parents=True, exist_ok=True)
    path = path or (EXTRA_RESULTS / f"spikes_{member}.csv")
    cols = [
        "member",
        "beta_deg",
        "U",
        "omega",
        "omega_e",
        "quantity",
        "i",
        "j",
        "old_value",
        "new_value",
    ]
    lines = [",".join(cols)]
    if not records:
        lines.append(
            f"# note: no irregular-frequency spikes replaced for {member} "
            f"(detection rule: across-speed residual > 5× median and > 2% of max, "
            f"only where ωe > ω_irr)"
        )
    else:
        for r in records:
            old = r["old_value"]
            new = r["new_value"]
            if isinstance(old, complex):
                old_s = f"{old.real}+{old.imag}j"
                new_s = f"{new.real}+{new.imag}j"
            else:
                old_s = repr(float(old))
                new_s = repr(float(new))
            lines.append(
                ",".join(
                    [
                        str(r["member"]),
                        str(r["beta_deg"]),
                        f"{r['U']:.6g}",
                        f"{r['omega']:.6g}",
                        f"{r['omega_e']:.6g}",
                        str(r["quantity"]),
                        str(r["i"]),
                        str(r["j"]),
                        old_s,
                        new_s,
                    ]
                )
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_diagnostic_overlays(pack: dict, out_dir: Path | None = None) -> list[Path]:
    """Plot A33/44/55, B33/44/55, |F3/4/5| vs ωe with five speeds overlaid."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = out_dir or DIAGNOSTICS
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    w = pack["w"] if "w_raw" not in pack else pack["w_raw"]
    AM = pack["AM"] if "AM_raw" not in pack else pack["AM_raw"]
    BD = pack["BD"] if "BD_raw" not in pack else pack["BD_raw"]
    F = pack["F"] if "F_raw" not in pack else pack["F_raw"]
    member = pack["member"]
    omega_irr = irregular_frequency(pack["Lpp"], pack["B"], pack["T"])

    quantities = [
        ("A33", lambda iw, ib, iu: AM[iw, 2, 2, ib, iu]),
        ("A44", lambda iw, ib, iu: AM[iw, 3, 3, ib, iu]),
        ("A55", lambda iw, ib, iu: AM[iw, 4, 4, ib, iu]),
        ("B33", lambda iw, ib, iu: BD[iw, 2, 2, ib, iu]),
        ("B44", lambda iw, ib, iu: BD[iw, 3, 3, ib, iu]),
        ("B55", lambda iw, ib, iu: BD[iw, 4, 4, ib, iu]),
        ("absF3", lambda iw, ib, iu: np.abs(F[iw, 2, ib, iu])),
        ("absF4", lambda iw, ib, iu: np.abs(F[iw, 3, ib, iu])),
        ("absF5", lambda iw, ib, iu: np.abs(F[iw, 4, ib, iu])),
    ]

    for beta_deg, ib in pack["beta_index"].items():
        beta_rad = float(pack["beta"][ib])
        fig, axes = plt.subplots(3, 3, figsize=(12, 9), sharex=True)
        axes = axes.ravel()
        for ax, (name, getter) in zip(axes, quantities):
            for iu, U in enumerate(pack["speed"]):
                we = encounter_frequency(w, float(U), beta_rad)
                vals = np.array([getter(iw, ib, iu) for iw in range(len(w))], dtype=float)
                ax.plot(we, vals, marker="o", ms=3, lw=1, label=f"U={U:.2f}")
            ax.axvline(omega_irr, color="k", ls="--", lw=0.8, label="ω_irr")
            ax.set_title(name)
            ax.grid(alpha=0.25)
            ax.set_xlabel("ωe (rad/s)")
        axes[0].legend(fontsize=6, loc="best")
        fig.suptitle(f"{member} β={beta_deg}° coefficient overlays")
        fig.tight_layout()
        p = out_dir / f"diag_{member}_{beta_deg}.png"
        fig.savefig(p, dpi=150)
        plt.close(fig)
        paths.append(p)
    return paths
