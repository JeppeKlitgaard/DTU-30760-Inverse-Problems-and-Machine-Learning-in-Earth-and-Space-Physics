import numba
import numpy as np
import numpy.typing as npt
from typing import Protocol, Any
from tqdm.auto import trange

EPS = 1e-14  # Small value to prevent division by zero


@numba.njit
def huber_weights(residuals: npt.NDArray, c: float) -> npt.NDArray:
    denom = np.maximum(np.abs(residuals), EPS)
    return np.minimum(c / denom, np.ones_like(residuals))


def median_absolute_deviation_scale_estimate(x: npt.NDArray) -> float:
    CORRECTION_FACTOR = 1 / 0.675  # = 1/(Φ^(-1)(0.75))
    median = np.median(x)
    mad = np.median(np.abs(x - median))

    return max(CORRECTION_FACTOR * mad, EPS)


class MEstimator(Protocol):
    def __call__(self, scaled_residuals: npt.NDArray, /) -> npt.NDArray: ...


class ScaleEstimator(Protocol):
    def __call__(self, residuals: npt.NDArray, /) -> float: ...


class HuberWeights:
    def __init__(self, c: float = 1.345):
        self.c = c  # Cut-off

    @staticmethod
    @numba.njit
    def _inner(residuals: npt.NDArray, c: float) -> npt.NDArray:
        denom = np.maximum(np.abs(residuals), EPS)
        return np.minimum(c / denom, np.ones_like(residuals))

    def __call__(self, residuals: npt.NDArray) -> npt.NDArray:
        return self._inner(residuals, c=self.c)


@numba.njit
def laplace_weights(residuals: npt.NDArray) -> npt.NDArray:
    denom = np.maximum(np.abs(residuals), EPS)
    return 1 / denom


@numba.njit
def _get_regularisation_weights(gamma: npt.NDArray, p: float, delta: float):
    ekblom = (gamma**2 + delta**2) ** (p / 2 - 1)
    return ekblom


@numba.njit
def _construct_A(
    G: npt.NDArray,
    w: npt.NDArray,
    H_k: tuple[npt.NDArray],
    alpha_k: tuple[float],
    delta_k: tuple[float],
    p_k: tuple[float],
    gamma_k: tuple[npt.NDArray],
):
    N_m = G.shape[1]
    K = len(H_k)

    reg_term = np.zeros((N_m, N_m))

    for k in numba.literal_unroll(range(K)):
        gamma = gamma_k[k]
        p = p_k[k]
        alpha = alpha_k[k]
        delta = delta_k[k]
        H = H_k[k]

        w_m = _get_regularisation_weights(gamma, p, delta)
        reg_term += (alpha**2) * H.T @ (w_m[:, None] * H)

    GTWG = G.T @ (w[:, None] * G)
    A = GTWG + reg_term

    return A


def irls_lp_multi_robust(
    G: npt.NDArray,  # Model
    d: npt.NDArray,  # Data
    H_k: tuple[npt.NDArray],  # Regularisers
    alpha_k: tuple[float],  # Reg. Parameters
    delta_k: tuple[float],  # Ekblom perturbations
    p_k: tuple[float],  # Norm orders
    robust_weighter: MEstimator = HuberWeights(),  # Robust weight function
    scale_estimator: ScaleEstimator = median_absolute_deviation_scale_estimate,  # Scale estimator
    N_iter: int = 1000,  # Maximum number of iterations
    diff_threshold: float = 1e-6,  # Tolerance for early stopping
    relative_diagonal_shift: float
    | None = None,  # Diagonal shift to improve conditioning
    return_iterates: bool = False,
    progress: bool = False,
    verbose: bool | int = False,
    freeze_scaler_after: int
    | None = None,  # Optionally freeze scale estimator after this many iterations
    check_condition: bool = False,  # Compute condition number and report it in output
    check_ekblom: bool = False,  # Optionally check and report Ekblom measures for debugging
) -> tuple[
    npt.NDArray,
    npt.NDArray,
    tuple[npt.NDArray],
    list[dict[str, Any]] | None,
    dict[str, Any],
]:
    """
    Solves a robust IRLS problem with an arbitrary number of regularisation terms
    and norm orders.

    Arguments:
        G:
            Forward/design matrix.
        d:
            Data vector (reshaped to 1D internally).
        H_k:
            Tuple of regularisation operators. Each operator maps model parameters
            to a regularisation domain for one term.
        alpha_k:
            Regularisation strengths for each term in H_k.
        delta_k:
            Ekblom perturbation constants for each regularisation term
            (stabilises LP-like weights).
        p_k:
            Norm orders for each regularisation term (for example 1.0 for L1-like behaviour).
        robust_weighter:
            Callable that maps scaled residuals to data weights.
        scale_estimator:
            Callable that estimates residual scale before robust weighting.
        N_iter:
            Maximum number of IRLS iterations.
        diff_threshold:
            Relative convergence threshold on model updates:
            ||m_new - m_old|| / (||m_new|| + eps).
        relative_diagonal_shift:
            Relative diagonal loading added to the system matrix each iteration:
            shift = relative_diagonal_shift * max(diag(A)).
            Helps numerical stability for ill-conditioned systems.
            Not used if None, otherwise it should be a small value like 1e-10 or smaller.
        return_iterates:
            If True, collect per-iteration diagnostics in iter_data.
        progress:
            If True, show a progress bar.
        verbose:
            If True, print per-iteration diagnostics.
            If it is an integer, only verbose every `i % verbose` iterations.
        freeze_scaler_after:
            If provided, update residual scale only for iterations i < freeze_scaler_after.
            After that, the most recent scale value is reused.
        check_condition:
            If True, compute and optionally warn about the condition number
            of the system matrix.
            Is expensive to calculate, so only use for debugging or if you suspect conditioning issues.
        check_ekblom:
            If True, compute and optionally report Ekblom measures for debugging.
    Returns:
        m: Estimated model vector.
        w: Final robust data weights.
        gamma_k: Final Ekblom measures.
        iter_data:
            Per-iteration diagnostics if return_iterates is True, otherwise None.
        info: Summary dictionary.
    """

    ### Initialisation
    N_d = G.shape[0]  # Number of data points
    N_m = G.shape[1]  # Number of model parameters
    K = len(H_k)  # Number of regularisation terms
    N_k = tuple(H.shape[0] for H in H_k)  # Sizes of regularisation domains

    ### Validation
    d_vec = d.reshape(-1)
    assert N_iter > 0
    if freeze_scaler_after is not None:
        assert freeze_scaler_after > 0
    assert d_vec.shape == (N_d,)

    assert len(alpha_k) == K
    assert len(delta_k) == K
    assert len(p_k) == K

    ### Initialisation
    m = np.zeros((N_m,))  # Model parameters
    w = np.ones((N_d,))  # Data weights
    gamma_k: tuple[npt.NDArray] = tuple(
        np.zeros((N_k[k],)) for k in range(K)
    )  # Ekblom measures

    if return_iterates:
        iter_data = []
    else:
        iter_data = None

    if progress:
        iter_ = trange(N_iter)
    else:
        iter_ = range(N_iter)

    threshold_check = np.inf  # Initialize threshold check for convergence
    residuals = np.zeros_like(np.inf)  # Initialize residuals
    scale = 1.0

    for i in iter_:
        ### Setup system
        A = _construct_A(G, w, H_k, alpha_k, delta_k, p_k, gamma_k)
        b = G.T @ (w * d_vec)

        ### Optionally shift diagonal to improve conditioning prior to solve
        if relative_diagonal_shift is not None:
            diag_shift = relative_diagonal_shift * np.max(np.diag(A))
            A[np.diag_indices_from(A)] += diag_shift

        ### Optional: Check condition number
        if check_condition:
            cond_number = np.linalg.cond(A)
            if cond_number > 1e12:
                print(
                    f"Warning: System is ill-conditioned at iteration {i} (cond={cond_number:.2e}). Results may be unstable."
                )

        ### Solve system
        # System is SPD, can use Cholesky for speedup
        m_new = scipy.linalg.solve(
            A, b, assume_a="pos"
        )  # Indicate that A is symmetric positive definite for speedup
        # m_new = scipy.linalg.solve(A, b)
        # m_new = np.linalg.solve(
        #     A, b,
        # )

        if not np.all(np.isfinite(m_new)):
            print(
                f"Non-finite values encountered in model update at iteration {i}. Terminating."
            )
            break

        ### Compute residuals and check convergence
        m_diff = m - m_new
        m_diff_L2 = np.linalg.norm(m_diff)

        m = m_new

        residuals = d_vec - G @ m

        ### Update robust weights
        if freeze_scaler_after is None or (
            freeze_scaler_after is not None and i < freeze_scaler_after
        ):
            scale = scale_estimator(residuals)

        w = robust_weighter(residuals / scale)

        ekblom_checks = []
        for k, H in enumerate(H_k):
            gamma_k[k][:] = H @ m

            if check_ekblom:
                ekblom_check = np.linalg.norm(gamma_k[k]) ** 2 / N_k[k]
                ekblom_checks.append(ekblom_check)

                if delta_k[k] ** 2 > (ekblom_check * 1e-3):
                    print(
                        f"Warning: Ekblom measure for regularisation term {k} is very small (||γ||^2/N={ekblom_check:.2e}), which may lead to numerical instability. Consider increasing delta_k[{k}]={delta_k[k]:.2e}."
                    )

        threshold_check = m_diff_L2 / (np.linalg.norm(m) + EPS)

        ### Debugging information
        if return_iterates:
            entry = {
                "i": i,
                "scale": scale,
                "threshold_check": threshold_check,
            }

            if check_condition:
                entry["cond"] = cond_number

            if check_ekblom:
                entry["ekblom_checks"] = ekblom_checks

            iter_data.append(entry)

        if verbose:
            if isinstance(verbose, int):
                if i % verbose != 0:
                    continue

            debug_msgs = [
                f"Iter {i}",
                f"m_diff_L2/||m||: {threshold_check:.2e}",
                f"max resid: {np.max(np.abs(residuals)):.2e}",
            ]
            if check_condition:
                debug_msgs.append(f"Condition number: {cond_number:.2e}")

            if check_ekblom:
                for k in range(K):
                    debug_msgs.append(
                        f"sqrt(||γ_{k}||_2^2 / N_k): {np.sqrt(ekblom_checks[k]):.2e}"
                    )
            print(", ".join(debug_msgs))

        ### Terminate if converged
        if threshold_check < diff_threshold:
            break

    info = {
        "iterations": i + 1,
        "relative_norm": threshold_check,
        "residuals": residuals,
    }

    return m, w, gamma_k, iter_data, info
