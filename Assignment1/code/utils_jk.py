import math

def format_number(x: float, sigfigs: int = 2) -> str:
    """
    Convert a number to MathText scientific notation.
    Example: 0.0053 -> r"5.3×10^{-3}"
    """
    if math.isnan(x):
        return r"$\mathrm{NaN}$"
    if math.isinf(x):
        return r"$\infty$" if x > 0 else r"$-\infty$"
    if x == 0:
        return r"$0$"

    exponent = int(math.floor(math.log10(abs(x))))
    mantissa = x / (10 ** exponent)

    # Keep requested significant figures in mantissa
    decimals = max(sigfigs - 1, 0)
    mantissa_str = f"{mantissa:.{decimals}f}".rstrip("0").rstrip(".")

    return rf"{mantissa_str}×10^{{{exponent}}}"
