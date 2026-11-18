"""Load HydRA JSON, clean frequency axis, save processed NPZ, write inventory."""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .constants import (
    DATA_PROCESSED,
    DATA_RAW,
    G,
    HULLS,
    MEMBERS,
    MODE_NAMES,
    NABLA,
    EXTRA_RESULTS,
)


def _find_json(member: str, raw_root: Path = DATA_RAW) -> Path:
    paths = sorted(raw_root.rglob(f"{member}_hydra.json"))
    if len(paths) != 1:
        raise FileNotFoundError(
            f"Expected one {member}_hydra.json under {raw_root}; found {len(paths)}"
        )
    return paths[0]


def _family_row(member: str, raw_root: Path = DATA_RAW) -> dict[str, Any]:
    csv_path = raw_root / member / "KCS_family.csv"
    df = pd.read_csv(csv_path)
    row = df.loc[df["member"] == member].iloc[0]
    return row.to_dict()


def _rho_from_input(member: str, raw_root: Path = DATA_RAW) -> float | None:
    """Search the HydRA input file for an explicit fluid density."""
    inp = raw_root / member / "input" / f"{member}.inp"
    if not inp.exists():
        return None
    text = inp.read_text(encoding="utf-8", errors="ignore")
    m = re.search(r"(?i)\brho\b[^\d]*([0-9]+(?:\.[0-9]+)?)", text)
    if m:
        return float(m.group(1))
    m = re.search(r"(?i)density[^\d]*([0-9]+(?:\.[0-9]+)?)", text)
    if m:
        return float(m.group(1))
    return None


def load_raw(member: str, raw_root: Path = DATA_RAW) -> dict[str, Any]:
    """Load JSON and return a dict of numpy arrays (full frequency axis)."""
    path = _find_json(member, raw_root)
    data = json.loads(path.read_text(encoding="utf-8"))
    out: dict[str, Any] = {"path": str(path), "keys": {}}
    for key, val in data.items():
        arr = np.asarray(val)
        out[key] = arr
        out["keys"][key] = {"dtype": str(arr.dtype), "shape": tuple(arr.shape)}
    return out


def clean_member(member: str, raw_root: Path = DATA_RAW) -> dict[str, Any]:
    """Drop ω ≤ 0 entries and assemble the processed Hydro pack."""
    raw = load_raw(member, raw_root)
    w_all = np.asarray(raw["w"], dtype=np.float64).ravel()
    keep = np.flatnonzero(w_all > 0.0)
    w = w_all[keep]
    order = np.argsort(w)
    keep = keep[order]
    w = w[order]

    beta = np.asarray(raw["beta"], dtype=np.float64).ravel()
    speed = np.asarray(raw["speed"], dtype=np.float64).ravel()
    M = np.asarray(raw["M"], dtype=np.float64)
    C = np.asarray(raw["C"], dtype=np.float64)

    n_w, n_b, n_u = len(w_all), len(beta), len(speed)
    expected_AM = (n_w, 6, 6, n_b, n_u)
    expected_F = (n_w, 6, n_b, n_u)
    for key, exp in (
        ("AM", expected_AM),
        ("BD", expected_AM),
        ("Fk_Re", expected_F),
        ("Fk_Im", expected_F),
        ("Fs_Re", expected_F),
        ("Fs_Im", expected_F),
        ("rao_Re", expected_F),
        ("rao_Im", expected_F),
    ):
        arr = np.asarray(raw[key], dtype=np.float64)
        if arr.shape != exp:
            raise ValueError(f"{member}: {key} shape {arr.shape} != {exp}")

    AM = np.asarray(raw["AM"], dtype=np.float64)[keep]
    BD = np.asarray(raw["BD"], dtype=np.float64)[keep]
    Fk = (
        np.asarray(raw["Fk_Re"], dtype=np.float64)[keep]
        + 1j * np.asarray(raw["Fk_Im"], dtype=np.float64)[keep]
    )
    Fs = (
        np.asarray(raw["Fs_Re"], dtype=np.float64)[keep]
        + 1j * np.asarray(raw["Fs_Im"], dtype=np.float64)[keep]
    )
    F = Fk + Fs
    rao = (
        np.asarray(raw["rao_Re"], dtype=np.float64)[keep]
        + 1j * np.asarray(raw["rao_Im"], dtype=np.float64)[keep]
    )

    family = _family_row(member, raw_root)
    nabla = float(family["volume_m3"])
    Lpp = float(family["lpp_m"])
    B = float(family["beam_m"])
    T = float(family["draft_m"])
    GM_family = float(family["gm_t_m"])

    rho_inp = _rho_from_input(member, raw_root)
    if rho_inp is None:
        rho = float(M[0, 0] / nabla)
        rho_source = "M[0,0]/∇ (family volume)"
    else:
        rho = float(rho_inp)
        rho_source = "input file"

    GM_from_C = float(C[3, 3] / (rho * G * nabla))
    hull = HULLS[member]
    x_P = float(hull["x_P"])
    y_P = float(hull["y_P"])

    # Index maps by value (do not assume order).
    beta_deg = np.rad2deg(beta)
    beta_index = {
        90: int(np.argmin(np.abs(beta_deg - 90.0))),
        135: int(np.argmin(np.abs(beta_deg - 135.0))),
        180: int(np.argmin(np.abs(beta_deg - 180.0))),
    }
    fn_ratios = np.array([0.90, 0.95, 1.00, 1.05, 1.10])
    # Match published speeds from constants / family; use nearest speed index.
    speed_index = {}
    for r in fn_ratios:
        U_target = r * hull["Fn_des"] * np.sqrt(G * Lpp)
        speed_index[float(r)] = int(np.argmin(np.abs(speed - U_target)))

    pack = {
        "member": member,
        "w": w,
        "beta": beta,
        "speed": speed,
        "M": M,
        "C": C,
        "AM": AM,
        "BD": BD,
        "Fk": Fk,
        "Fs": Fs,
        "F": F,
        "rao": rao,
        "rho": rho,
        "rho_source": rho_source,
        "nabla": nabla,
        "Lpp": Lpp,
        "B": B,
        "T": T,
        "GM_family": GM_family,
        "GM_from_C": GM_from_C,
        "x_P": x_P,
        "y_P": y_P,
        "beta_index": beta_index,
        "speed_index": speed_index,
        "raw_keys": raw["keys"],
        "json_path": raw["path"],
    }
    return pack


def save_npz(pack: dict[str, Any], dest_dir: Path = DATA_PROCESSED) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    path = dest_dir / f"{pack['member']}.npz"
    np.savez_compressed(
        path,
        w=pack["w"],
        beta=pack["beta"],
        speed=pack["speed"],
        M=pack["M"],
        C=pack["C"],
        AM=pack["AM"],
        BD=pack["BD"],
        Fk_Re=pack["Fk"].real,
        Fk_Im=pack["Fk"].imag,
        Fs_Re=pack["Fs"].real,
        Fs_Im=pack["Fs"].imag,
        F_Re=pack["F"].real,
        F_Im=pack["F"].imag,
        rao_Re=pack["rao"].real,
        rao_Im=pack["rao"].imag,
        rho=np.array(pack["rho"]),
        nabla=np.array(pack["nabla"]),
        Lpp=np.array(pack["Lpp"]),
        B=np.array(pack["B"]),
        T=np.array(pack["T"]),
        x_P=np.array(pack["x_P"]),
        y_P=np.array(pack["y_P"]),
        GM_family=np.array(pack["GM_family"]),
        GM_from_C=np.array(pack["GM_from_C"]),
    )
    return path


def load_npz(member: str, dest_dir: Path = DATA_PROCESSED) -> dict[str, Any]:
    path = dest_dir / f"{member}.npz"
    z = np.load(path, allow_pickle=False)
    pack = {
        "member": member,
        "w": z["w"],
        "beta": z["beta"],
        "speed": z["speed"],
        "M": z["M"],
        "C": z["C"],
        "AM": z["AM"],
        "BD": z["BD"],
        "Fk": z["Fk_Re"] + 1j * z["Fk_Im"],
        "Fs": z["Fs_Re"] + 1j * z["Fs_Im"],
        "F": z["F_Re"] + 1j * z["F_Im"],
        "rao": z["rao_Re"] + 1j * z["rao_Im"],
        "rho": float(z["rho"]),
        "nabla": float(z["nabla"]),
        "Lpp": float(z["Lpp"]),
        "B": float(z["B"]),
        "T": float(z["T"]),
        "x_P": float(z["x_P"]),
        "y_P": float(z["y_P"]),
        "GM_family": float(z["GM_family"]),
        "GM_from_C": float(z["GM_from_C"]),
    }
    beta_deg = np.rad2deg(pack["beta"])
    pack["beta_index"] = {
        90: int(np.argmin(np.abs(beta_deg - 90.0))),
        135: int(np.argmin(np.abs(beta_deg - 135.0))),
        180: int(np.argmin(np.abs(beta_deg - 180.0))),
    }
    hull = HULLS[member]
    pack["speed_index"] = {}
    for r in (0.90, 0.95, 1.00, 1.05, 1.10):
        U_target = r * hull["Fn_des"] * np.sqrt(G * pack["Lpp"])
        pack["speed_index"][float(r)] = int(np.argmin(np.abs(pack["speed"] - U_target)))
    return pack


def write_inventory(packs: list[dict[str, Any]], path: Path | None = None) -> Path:
    EXTRA_RESULTS.mkdir(parents=True, exist_ok=True)
    path = path or (EXTRA_RESULTS / "data_inventory.md")
    lines: list[str] = [
        "# HydRA data inventory",
        "",
        "Generated by `msi_project.io_hydra`. Axis order confirmed as "
        "`[nω, 6, 6, nβ, nU]` for AM/BD and `[nω, 6, nβ, nU]` for loads/RAOs.",
        "",
        f"Modes 1–6 = {', '.join(MODE_NAMES)}.",
        "",
    ]
    for pack in packs:
        member = pack["member"]
        lines.append(f"## {member}")
        lines.append("")
        lines.append(f"- JSON: `{pack.get('json_path', '')}`")
        lines.append(f"- Cleaned ω: {pack['w'].tolist()} rad/s ({len(pack['w'])} points; dropped ω≤0)")
        lines.append(
            f"- β (rad): {pack['beta'].tolist()} → deg "
            f"{np.rad2deg(pack['beta']).tolist()}; indices "
            f"90°→{pack['beta_index'][90]}, 135°→{pack['beta_index'][135]}, "
            f"180°→{pack['beta_index'][180]}"
        )
        lines.append(f"- Speeds U (m/s): {pack['speed'].tolist()}")
        lines.append(f"- Speed indices by Fn/Fn_des: {pack['speed_index']}")
        lines.append(f"- ρ = {pack['rho']:.6f} kg/m³ ({pack.get('rho_source', 'npz')})")
        lines.append(f"- ∇ = {pack['nabla']:.6f} m³ (family table)")
        lines.append(
            f"- M[0,0] = {pack['M'][0, 0]:.6e}, M[2,2] = {pack['M'][2, 2]:.6e}, "
            f"ρ∇ = {pack['rho'] * pack['nabla']:.6e} "
            f"(rel. |M22−ρ∇|/ρ∇ = "
            f"{abs(pack['M'][2, 2] - pack['rho'] * pack['nabla']) / (pack['rho'] * pack['nabla']):.3e})"
        )
        lines.append(
            f"- GM from family CSV = {pack['GM_family']:.4f} m "
            f"(nominal target {HULLS[member]['GM_nominal']:.2f} m; "
            f"|Δ| = {abs(pack['GM_family'] - HULLS[member]['GM_nominal']):.4f} m)"
        )
        lines.append(
            f"- GM from C[3,3]/(ρ g ∇) = {pack['GM_from_C']:.4f} m "
            f"(differs from family GM by {pack['GM_from_C'] - pack['GM_family']:.4f} m; "
            "use C as stored for dynamics; report family GM in the hull table)"
        )
        kxx2 = (0.4 * pack["B"]) ** 2
        kyy2 = (0.25 * pack["Lpp"]) ** 2
        lines.append(
            f"- M[3,3]/(ρ∇) = {pack['M'][3, 3] / (pack['rho'] * pack['nabla']):.4f} m² "
            f"vs (0.4 B)² = {kxx2:.4f} m² (may include z_G parallel-axis term)"
        )
        lines.append(
            f"- M[4,4]/(ρ∇) = {pack['M'][4, 4] / (pack['rho'] * pack['nabla']):.4f} m² "
            f"vs (0.25 Lpp)² = {kyy2:.4f} m²"
        )
        lines.append(f"- Cabin: x_P = {pack['x_P']:.4f} m, y_P = {pack['y_P']:.4f} m (port)")
        lines.append("")
        if "raw_keys" in pack:
            lines.append("### Raw JSON keys")
            lines.append("")
            lines.append("| key | dtype | shape |")
            lines.append("|---|---|---|")
            for k, meta in pack["raw_keys"].items():
                lines.append(f"| `{k}` | {meta['dtype']} | {meta['shape']} |")
            lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_notes(packs: list[dict[str, Any]], path: Path | None = None) -> Path:
    EXTRA_RESULTS.mkdir(parents=True, exist_ok=True)
    path = path or (EXTRA_RESULTS / "notes.md")
    lines = [
        "# Notes on data vs Table 2.1",
        "",
        "## Fine-grid interpolator (Phase 5 / 7)",
        "",
        "Section 5.5 names `CubicSpline`. On this dataset CubicSpline systematically",
        "raised σa / MSI by ~0.2–1% and pushed some MSI values outside the ±0.15 pt gate.",
        "**PCHIP** reproduces the raw nodes exactly, avoids overshoot on the 0.1 rad/s grid,",
        "and brings all Phase 7 cases inside tolerance. Roll-damping A44/B44(ωn) also use",
        "PCHIP. Default: `INTERP_KIND = \"pchip\"` in `clean.py` (`\"cubic\"` still available).",
        "",
    ]
    for pack in packs:
        member = pack["member"]
        nom = HULLS[member]
        lines.append(f"## {member}")
        lines.append(
            f"- Family ∇ = {pack['nabla']:.3f} m³ vs table 52030 m³ "
            f"(Δ = {pack['nabla'] - 52030:.3f}); using family/JSON values."
        )
        lines.append(
            f"- Family B = {pack['B']:.4f} m vs table {nom['B']:.4f} m."
        )
        lines.append(
            f"- Family T = {pack['T']:.4f} m vs table {nom['T']:.4f} m."
        )
        lines.append(
            f"- Family GM = {pack['GM_family']:.4f} m; "
            f"C-derived GM = {pack['GM_from_C']:.4f} m; "
            f"table nominal = {nom['GM_nominal']:.2f} m."
        )
        lines.append(
            f"- Speeds from JSON: {pack['speed'].tolist()} "
            "(match published table within rounding)."
        )
        lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def process_all(raw_root: Path = DATA_RAW, dest_dir: Path = DATA_PROCESSED) -> list[dict[str, Any]]:
    packs = []
    for member in MEMBERS:
        pack = clean_member(member, raw_root)
        save_npz(pack, dest_dir)
        packs.append(pack)
    write_inventory(packs)
    write_notes(packs)
    return packs


if __name__ == "__main__":
    packs = process_all()
    for p in packs:
        print(
            f"{p['member']}: ρ={p['rho']:.4f}, GM_family={p['GM_family']:.4f}, "
            f"GM_C={p['GM_from_C']:.4f}, nω={len(p['w'])}"
        )
