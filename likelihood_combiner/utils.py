import numpy as np

__all__ = [
    'get_sigmav_range',
    'round_sigmav_range',
    'progress_bar'
]

def get_sigmav_range(sigmav_min=1e-28, sigmav_max=1e-18, sigmav_nPoints=1001, precision=3):

    # Constructing the sigmav range and spacing.
    sigmav_min = -(np.abs(np.floor(np.log10(np.abs(float(sigmav_min)))).astype(int))).astype(int)
    sigmav_max = -(np.abs(np.floor(np.log10(np.abs(float(sigmav_max)))).astype(int))).astype(int)
    sigmav_nPoints = int(sigmav_nPoints)
    sigmav_range = np.logspace(sigmav_min, sigmav_max, sigmav_nPoints, endpoint=True, dtype=np.float32)
    return round_sigmav_range(sigmav_range, precision)

def round_sigmav_range(sigmav_range, precision=3):
    # Round of the third (default) digit to avoid an interpolation for the GloryDuck files.
    exponent = (np.abs(np.floor(np.log10(np.abs(sigmav_range))).astype(int))+precision).astype(int)
    for i,e in enumerate(exponent):
        sigmav_range[i] = np.around(sigmav_range[i],decimals=e)
    return sigmav_range

def progress_bar(current_value, total):
    increments = 50
    percentual = ((current_value/ total) * 100)
    i = np.int(percentual // (100 / increments ))
    percentual = np.int(percentual)
    text = "\r[{0: <{1}}] {2}%".format('=' * i, increments, percentual)
    print(text, end="\n" if percentual == 100 else "")

