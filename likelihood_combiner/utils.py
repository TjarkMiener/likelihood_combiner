"""
utils.py
========
Helper functions, which might be useful.
"""

import numpy as np

__all__ = [
    'get_sigmav_range',
    'round_sigmav_range',
    'progress_bar'
]

def get_sigmav_range(sigmav_min=1e-28, sigmav_max=1e-18, sigmav_nPoints=1001, precision=3):
    """
    Retrieve the sigmav range in log-spacing.

    Parameters
    ----------
    sigmav_min: float
        minimum value of the sigmav range.
    sigmav_max: float
        maximum value of the sigmav range.
    sigmav_nPoints: int
        number of points in the sigmav range.
    precision: int
        precision of the returning sigmav range.

    Returns
    -------
    sigmav_range: `numpy.ndarray of type numpy.float32`
        sigmav range (ascending).
    """

    # Constructing the sigmav range with log-spacing.
    sigmav_min = -(np.abs(np.floor(np.log10(np.abs(float(sigmav_min)))).astype(int))).astype(int)
    sigmav_max = -(np.abs(np.floor(np.log10(np.abs(float(sigmav_max)))).astype(int))).astype(int)
    sigmav_nPoints = int(sigmav_nPoints)
    sigmav_range = np.logspace(sigmav_min, sigmav_max, sigmav_nPoints, endpoint=True, dtype=np.float32)
    return round_sigmav_range(sigmav_range, precision)

def round_sigmav_range(sigmav_range, precision=3):
    """
    Round the sigmav range.

    Parameters
    ----------
    sigmav_range: `numpy.ndarray of type numpy.float32`
        sigmav range (ascending).
    precision: int
        precision of the returning sigmav range.

    Returns
    -------
    sigmav_range: `numpy.ndarray of type numpy.float32`
        rounded sigmav range (ascending).
    """

    # Round of the third (default) digit to avoid an interpolation for the GloryDuck files.
    exponent = (np.abs(np.floor(np.log10(np.abs(sigmav_range))).astype(int))+precision).astype(int)
    for i,e in enumerate(exponent):
        sigmav_range[i] = np.around(sigmav_range[i],decimals=e)
    return sigmav_range

def progress_bar(current_value, total):
    """
    Print progress bar in parallel processing mode.

    Parameters
    ----------
    current_value: int
        current status of the process.
    total: int
        total number of simulation.
    """

    increments = 50
    percentual = ((current_value/ total) * 100)
    i = np.int(percentual // (100 / increments ))
    percentual = np.int(percentual)
    text = "\r[{0: <{1}}] {2}%".format('=' * i, increments, percentual)
    print(text, end="\n" if percentual == 100 else "")

