"""O'Hanlon & McCauley Motion Sickness Incidence model."""
from __future__ import annotations

import numpy as np
from scipy.stats import norm

from .constants import G


def mu_msi(we_bar: float) -> float:
    """Frequency factor μ_MSI from mean encounter frequency (rad/s)."""
    return float(-0.819 + 2.32 * (np.log10(we_bar)) ** 2)


def msi_percent(abar: float, we_bar: float, g: float = G) -> float:
    """MSI in percent of people who vomit within 2 hours.

    Φ is the standard normal CDF (scipy.stats.norm.cdf), NOT erf alone.
    """
    mu = mu_msi(we_bar)
    z = (np.log10(abar / g) - mu) / 0.4
    return float(100.0 * norm.cdf(z))


def from_sigma(sigma_a: float, we_bar: float, g: float = G) -> dict[str, float]:
    """Convert acceleration std-dev and mean encounter frequency to MSI outputs."""
    abar = float(np.sqrt(2.0 / np.pi) * sigma_a)
    mu = mu_msi(we_bar)
    msi = msi_percent(abar, we_bar, g=g)
    return {
        "sigma_a": float(sigma_a),
        "abar": abar,
        "we_bar": float(we_bar),
        "mu_msi": mu,
        "msi_pct": msi,
    }
