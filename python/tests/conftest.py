"""Shared helpers for the test suite."""
import numpy as np


def grid(res, pole_points=False):
    """Regular global grid; with ``pole_points`` the rows sit on -90..90 (NCEP style)."""
    if pole_points:
        lat = np.arange(-90, 90 + res / 2, res)
    else:
        lat = np.arange(-90 + res / 2, 90, res)
    lon = np.arange(0, 360, res) + (0 if pole_points else res / 2)
    return lat, lon
