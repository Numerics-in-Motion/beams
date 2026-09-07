"""Utility helpers: config loading, JSON encoding, logging."""
from __future__ import annotations

import json
import os
import time
from contextlib import contextmanager

import numpy as np


def load_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


class NumpyJSONEncoder(json.JSONEncoder):
    """JSON encoder that understands numpy scalars / arrays."""

    def default(self, obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def save_json(obj, path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False, cls=NumpyJSONEncoder)


@contextmanager
def timer():
    """Context manager returning elapsed wall time in seconds via .value."""

    class _T:
        value = 0.0

    t = _T()
    start = time.perf_counter()
    try:
        yield t
    finally:
        t.value = time.perf_counter() - start


def area_from_effective_depth(section: dict) -> float:
    """A nominal reference area implied by the configured section depth.

    Only used as a fallback; the real per-structure area is set by the
    material-volume matching in geometry.py.
    """
    d = float(section["effective_depth_m"])
    if section.get("type", "circular") == "circular":
        return np.pi * (d / 2.0) ** 2
    # square section fallback
    return d * d


def circular_section_properties(area: float) -> tuple[float, float, float]:
    """Return (radius, moment_of_inertia, outer_fibre_c) for a solid circular
    cross-section of the given area.

    A = pi r^2  ->  r = sqrt(A/pi)
    I = pi r^4 / 4 = A^2 / (4 pi)
    c = r (distance from neutral axis to outer fibre)
    """
    r = np.sqrt(area / np.pi)
    inertia = area ** 2 / (4.0 * np.pi)
    return r, inertia, r
