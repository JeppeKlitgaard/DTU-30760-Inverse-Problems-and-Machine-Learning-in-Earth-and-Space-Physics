import numpy as np


def l_curve_corner_aster(rho, eta, reg_param):
    # Parameter Estimation and Inverse Problems, 2nd edition, 2011
    # by R. Aster, B. Borchers, C. Thurber
    #
    # reg_corner, ireg_corner, kappa = l_curve_corner(rho,eta,reg_param)
    # returns l curve corner estimated using a maximum curvature (kappa) estimation
    # in log-log space
    # rho is the misfit and eta is the model norm or seminorm
    #
    # INPUT
    #   rho       - misfit
    #   eta       - model norm or seminorm
    #   reg_param - the regularization parameter
    #
    # OUTPUT
    #   reg_corner  - the value of reg_param with maximum curvature
    #   ireg_corner - the index of the value in reg_param with maximum curvature
    #   kappa       - the curvature for each reg_param
    #
    # transform rho and eta into log-log space
    rho = np.array(rho, dtype=np.float64)
    eta = np.array(eta, dtype=np.float64)
    reg_param = np.array(reg_param, dtype=np.float64)

    x = np.log10(rho)
    y = np.log10(eta)

    # Triangular/circumscribed circle simple approximation to curvature
    # (after Roger Stafford)

    # the series of points used for the triangle/circle
    x1 = x[0:-2]
    x2 = x[1:-1]
    x3 = x[2:]
    y1 = y[0:-2]
    y2 = y[1:-1]
    y3 = y[2:]

    # the side lengths for each triangle
    a = np.sqrt((x3 - x2) ** 2 + (y3 - y2) ** 2)
    b = np.sqrt((x1 - x3) ** 2 + (y1 - y3) ** 2)
    c = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

    s = (a + b + c) / 2  # semi-perimeter

    # the radius of each circle
    sa = s - a
    sb = s - b
    sc = s - c

    C_div = 4 * np.sqrt((s * (sa) * (sb) * (sc)))
    R = (a * b * c) / C_div

    # The curvature for each estimate for each value which is
    # the reciprocal of its circumscribed radius. Since there aren't circles for
    # the end points they have no curvature
    kappa = np.hstack((0, 1 / R, 0)).T
    ireg_corner = np.nanargmax(abs(kappa))
    reg_corner = reg_param[ireg_corner]

    return reg_corner, ireg_corner, kappa


# LLM Generated, Gemini 3.1 Pro, 2026-03-10
def l_curve_corner_strided_menger(rho, eta, alphas, stride=2):
    """
    Finds the L-curve corner using normalized, strided Menger curvature.

    Parameters
    ----------
    rho : array_like
        The residual norms (misfit).
    eta : array_like
        The model norms (seminorm).
    alphas : array_like
        The regularization parameters.
    stride : int
        The step size for forming triangles. stride=1 uses adjacent points.
        stride > 1 helps skip microscopic numerical jitter from solvers.

    Returns
    -------
    best_alpha : float
        The regularization parameter at the point of maximum curvature.
    best_idx : int
        The index of the optimal alpha.
    kappa : ndarray
        The computed curvature array (padded with zeros at the boundaries).
    """
    rho = np.array(rho, dtype=np.float64)
    eta = np.array(eta, dtype=np.float64)
    alphas = np.array(alphas, dtype=np.float64)
    N = len(alphas)

    if N <= 2 * stride:
        raise ValueError("Not enough points to compute curvature with the given stride.")

    ### Transform and normalise
    x = np.log10(rho)
    y = np.log10(eta)
    x_range = x.max() - x.min()
    y_range = y.max() - y.min()

    # Avoid division by zero in the degenerate case of a perfectly constant norm
    if x_range == 0: x_range = 1.0
    if y_range == 0: y_range = 1.0
    x_norm = (x - x.min()) / x_range
    y_norm = (y - y.min()) / y_range

    ### Stride
    # p1: points before, p2: center points (where curvature is evaluated), p3: points after
    x1 = x_norm[:-2 * stride]
    x2 = x_norm[stride:-stride]
    x3 = x_norm[2 * stride:]

    y1 = y_norm[:-2 * stride]
    y2 = y_norm[stride:-stride]
    y3 = y_norm[2 * stride:]

    ### Side lengths and area for Menger curvature
    a = np.sqrt((x3 - x2)**2 + (y3 - y2)**2)  # Distance p2 to p3
    b = np.sqrt((x1 - x3)**2 + (y1 - y3)**2)  # Distance p1 to p3
    c = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)  # Distance p1 to p2

    s = (a + b + c) / 2.0  # Area
    # np.maximum prevents invalid sqrt for near-degenerate floating point triangles
    area = np.sqrt(np.maximum(s * (s - a) * (s - b) * (s - c), 0.0))

    ### Menger curvature: κ = 4 * Area / (a * b * c)
    with np.errstate(divide='ignore', invalid='ignore'):
        kappa_inner = 4.0 * area / (a * b * c)
        # If points are collinear or coincident, a*b*c = 0. Set curvature to 0.
        kappa_inner = np.nan_to_num(kappa_inner, nan=0.0, posinf=0.0)

    # Pad to match the original array length and align indices correctly
    pad = np.zeros(stride)
    kappa = np.concatenate((pad, kappa_inner, pad))

    # Find the point of maximum curvature
    best_idx = np.argmax(kappa)
    best_alpha = alphas[best_idx]

    return best_alpha, best_idx, kappa


# LLM Generated, Gemini 3.1 Pro, 2026-03-10
def l_curve_corner_chord_distance(rho, eta, reg_param):
    """
    Finds the L-curve corner using the geometric chord-distance method.
    This method is highly robust to IRLS noise and aspect-ratio scaling.

    Parameters
    ----------
    rho : array_like
        The residual norms (misfit).
    eta : array_like
        The model norms (seminorm).
    reg_param : array_like
        The regularization parameters.

    Returns
    -------
    reg_corner : float
        The regularization parameter at the point furthest from the chord.
    ireg_corner : int
        The index of the optimal regularization parameter.
    distances : ndarray
        The computed chord distances for each reg_param.
    """
    rho = np.array(rho, dtype=np.float64)
    eta = np.array(eta, dtype=np.float64)
    reg_param = np.array(reg_param, dtype=np.float64)

    ### Transform and normalise
    x = np.log10(rho)
    y = np.log10(eta)
    x_range = x.max() - x.min()
    y_range = y.max() - y.min()

    # Avoid division by zero in the degenerate case of a perfectly constant norm
    if x_range == 0: x_range = 1.0
    if y_range == 0: y_range = 1.0

    x_norm = (x - x.min()) / x_range
    y_norm = (y - y.min()) / y_range

    ### Define the start and end points of the curve
    p_start_x = x_norm[0]
    p_start_y = y_norm[0]

    p_end_x = x_norm[-1]
    p_end_y = y_norm[-1]

    # Vector of the line connecting endpoints
    line_vec_x = p_end_x - p_start_x
    line_vec_y = p_end_y - p_start_y
    line_len = np.sqrt(line_vec_x**2 + line_vec_y**2)

    ### Calculate perpendicular distance from each point to the line (vectorized)
    # 2D cross product magnitude gives the area of parallelogram;
    # dividing by base length gives the height (perpendicular distance)
    cross_prod = (x_norm - p_start_x) * line_vec_y - (y_norm - p_start_y) * line_vec_x

    with np.errstate(divide='ignore', invalid='ignore'):
        distances = np.abs(cross_prod) / line_len
        distances = np.nan_to_num(distances, nan=0.0, posinf=0.0)

    # Find the point furthest from the chord
    ireg_corner = np.argmax(distances)
    reg_corner = reg_param[ireg_corner]

    return reg_corner, ireg_corner, distances
