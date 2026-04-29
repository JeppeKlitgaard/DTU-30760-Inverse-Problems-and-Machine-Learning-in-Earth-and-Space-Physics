# %% [markdown]
# # Imports

# %%
# JAX setup
import jax

# Enable 64-bit precision in JAX
jax.config.update("jax_enable_x64", True)

# Which accelerator to use
# jax.config.update("jax_platforms", "cpu")
jax.config.update("jax_platforms", "cuda,cpu")


# %%
import itertools
from typing import Literal

import cmcrameri.cm as cmc
import equinox as eqx
import gstools as gs
import jax.numpy as jnp
import joblib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from mpl_toolkits.axes_grid1 import ImageGrid, make_axes_locatable
from nanopinv._typing import Array, Float
from nanopinv.distribution import MultivariateNormalCholesky
from nanopinv.physics.eikonal import build_travel_time_points
from nanopinv.sampler import (
    ExtendedMetropolisChain,
    ParallelTemperingSampler,
    ProposalDistribution,
    initialize_betas,
)
from nanopinv.types import ObservationsUnivariate
from nanopinv.utils import StatefulRNGKey, make_pytree_spec
from nanopinv.variance import get_distance_matrix
from scipy.io import loadmat
from pathlib import Path

# %%
sns.set_theme(context="paper", style="whitegrid", rc={"figure.dpi": 150})
skey = StatefulRNGKey(0)
mem = joblib.Memory(".joblib_cache", verbose=0)

ASSIGNMENT_DIR = Path().resolve()
REPORT_DIR = ASSIGNMENT_DIR / "report"
EXPORT_DIR = REPORT_DIR / "export"
EXPORT_DIR.mkdir(exist_ok=True)
print(ASSIGNMENT_DIR)

A4_WIDTH_IN = 210 / 25.4
_AVAILABLE_WIDTH_MM_NORMAL = 160  # With steal-margin
_AVAILABLE_WIDTH_MM_EXTRA = 170  # With steal-margin
WIDTH_IN_NORMAL = _AVAILABLE_WIDTH_MM_NORMAL / 25.4
WIDTH_IN_EXTRA = _AVAILABLE_WIDTH_MM_EXTRA / 25.4
EXPORT_KWARGS = {
    "bbox_inches": "tight",
    "pad_inches": 0.025,
    # "dpi": 600,
    "dpi": 300,
}

# %% [markdown]
# # Load data

# %%
# Given by assignment
data_raw = loadmat("data/BHRS_2D.mat")
data_obs = data_raw["traveltimes"].squeeze()
data_std = data_raw["traveltimes_std"].squeeze()
ndata, ndim = data_raw["sources"].shape
sources_raw = data_raw["sources"].squeeze()
receivers_raw = data_raw["receivers"].squeeze()

# %%
# sources and receivers are 3D, but third dimension is all zeros
assert np.allclose(sources_raw[:, 2], 0.0)
assert np.allclose(receivers_raw[:, 2], 0.0)

# Remove third dimension to make it 2D problem, I guess this should have been done in the starting notebook...?
ndim = 2
sources = sources_raw[:, :2]
receivers = receivers_raw[:, :2]
ndata = sources.shape[0]
assert receivers.shape[0] == ndata

# %%
# Observations
data_obs = data_raw["traveltimes"].squeeze()
data_std = data_raw["traveltimes_std"].squeeze().astype(np.float64)

# We have univariate observations (according to the assumptions of our a priori model),
# so we can use optimized Observations class here
assert (data_std == 2.0).all()
obs = ObservationsUnivariate(data_obs=data_obs, data_std=2.0)

# %%
# Calculate extents for data grid
_data_grid_all_x = np.hstack((sources[:, 0], receivers[:, 0]))
_data_grid_all_y = np.hstack((sources[:, 1], receivers[:, 1]))

data_grid_min_x = np.min(_data_grid_all_x)
data_grid_max_x = np.max(_data_grid_all_x)
data_grid_min_y = np.min(_data_grid_all_y)
data_grid_max_y = np.max(_data_grid_all_y)

# %%
# Discretise and construct model grid
model_grid_dx = 0.25
model_grid_dy = model_grid_dx

model_grid_x = np.arange(
    data_grid_min_x, data_grid_max_x + model_grid_dx, model_grid_dx
)
model_grid_y = np.arange(
    data_grid_min_y, data_grid_max_y + model_grid_dy, model_grid_dy
)
model_grid_X, model_grid_Y = np.meshgrid(model_grid_x, model_grid_y, indexing="ij")
model_grid_shape = model_grid_X.shape  # (Nx, Ny)
model_grid_extent = [
    np.min(model_grid_x),
    np.max(model_grid_x),
    np.min(model_grid_y),
    np.max(model_grid_y),
]

# %%
# Print info for use in report
N_d = len(data_obs)
N_sources = len(np.unique(sources))
N_receivers = len(np.unique(receivers))

print(f"{N_d=}")
print(f"{N_sources=}")
print(f"{N_receivers=}")
print(f"y ∈ [{data_grid_min_y}, {data_grid_max_y}]")
print(f"x_ptp = {data_grid_max_x - data_grid_min_x}")
print(f"N_x = {model_grid_shape[0]}")
print(f"N_y = {model_grid_shape[1]}")
print(f"N_x × N_y = {model_grid_shape[0] * model_grid_shape[1]}")

# %%
# Load and preprocess prior data
data_empirical_prior_raw = loadmat("data/BHRS_prior.mat")

empirical_prior_model_grid_x = data_empirical_prior_raw["x"].squeeze()
empirical_prior_model_grid_y = data_empirical_prior_raw["y"].squeeze()
empirical_prior_model_grid_shape = (
    empirical_prior_model_grid_x.size,
    empirical_prior_model_grid_y.size,
)

empirical_prior_samples = [
    data_empirical_prior_raw["prior_0"],
    data_empirical_prior_raw["prior_1"],
    data_empirical_prior_raw["prior_2"],
]
empirical_prior_samples_arr = np.asarray(empirical_prior_samples)

empirical_prior_mean = np.mean(empirical_prior_samples)
print(f"{empirical_prior_mean=}")

# %%
# Plotting function for field
_PLOT_Y_TOP = 1.0


def _get_ax(ax=None, figsize=(4, 8)):
    """Creates a default axis if none is provided."""
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    return ax


def format_ax(ax, title=None, extent=None):
    """Applies common labels, grids, and precise title placement."""
    ax.set_xlabel("$x$ [m]")
    ax.set_ylabel("$y$ [m]")

    ax.set_aspect("equal")
    if extent is not None:
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])

    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.grid(True, alpha=0.5)
    if title is not None:
        ax.set_title(title, x=0.0, y=_PLOT_Y_TOP, size=12, horizontalalignment="left")


def add_legend(ax):
    """Standardized legend explicitly placed at top right."""
    ax.legend(
        loc="lower right",
        bbox_to_anchor=(1.0, _PLOT_Y_TOP),
        frameon=True,
        facecolor=plt.cm.binary(0.3),
    )


def plot_field(
    field,
    extent=model_grid_extent,
    ax=None,
    cmap=cmc.batlow,
    plot_cbar: bool = True,
    cbar_label: str = "Value",
    cbar_location: Literal["left", "right", "top", "bottom"] = "bottom",
    plot_grid_points: bool = True,
    grid_points=(model_grid_X, model_grid_Y),
    plot_water_table: bool = True,
    water_table_elevation: float = 844.0,
    auto_truncate: bool = True,
):
    ax = _get_ax(ax)

    # Imshow Field
    im = ax.imshow(field.T, origin="lower", extent=extent, cmap=cmap)

    if plot_cbar:
        # Dynamically attach the colorbar to the drawn axes
        divider = make_axes_locatable(ax)
        orientation = "horizontal" if cbar_location in ["top", "bottom"] else "vertical"
        cax = divider.append_axes(cbar_location, size="3%", pad=0.5)
        plt.colorbar(im, cax=cax, orientation=orientation, label=cbar_label)

    # Overlays
    if plot_grid_points and grid_points is not None:
        gx, gy = grid_points
        if auto_truncate:
            mask = gy < extent[3]
            gx, gy = gx[mask], gy[mask]
        ax.scatter(gx, gy, color="k", s=1, label="Model grid", alpha=0.5)

    if plot_water_table:
        ax.axhline(water_table_elevation, color="k", ls="--", label="Water Table")

    format_ax(ax, extent=extent)
    return ax


def plot_sender_receiver(
    sources=sources,
    receivers=receivers,
    extent=model_grid_extent,
    ax=None,
    auto_truncate: bool = True,
):
    ax = _get_ax(ax)

    if auto_truncate:
        s_mask = sources[:, 1] < extent[3]
        r_mask = receivers[:, 1] < extent[3]
        sources = sources[s_mask]
        receivers = receivers[r_mask]

    ax.scatter(
        receivers[:, 0],
        receivers[:, 1],
        color=cmc.buda(0.999),
        label="Receivers",
        marker="x",
        lw=0.75,
    )
    ax.scatter(
        sources[:, 0],
        sources[:, 1],
        color=cmc.buda(0),
        label="Sources",
        marker="o",
        facecolors="none",
        lw=0.75,
    )
    return ax


def plot_observations(
    sources=sources,
    receivers=receivers,
    observations=data_obs,
    extent=model_grid_extent,
    ax=None,
    cmap=cmc.buda,
    plot_cbar: bool = True,
    cbar_label: str = "Travel Time [ns]",
    cbar_location: Literal["left", "right", "top", "bottom"] = "bottom",
    auto_truncate: bool = True,
):
    ax = _get_ax(ax)

    if auto_truncate:
        joint_mask = (sources[:, 1] < extent[3]) & (receivers[:, 1] < extent[3])
        sources = sources[joint_mask]
        receivers = receivers[joint_mask]
        observations = observations[joint_mask]

    obs_min, obs_max = np.min(observations), np.max(observations)
    obs_range = obs_max - obs_min if obs_max > obs_min else 1.0

    for (sx, sy), (rx, ry), obs in zip(sources, receivers, observations):
        color_fraction = (obs - obs_min) / obs_range
        ax.plot([sx, rx], [sy, ry], color=cmap(color_fraction), alpha=0.5, lw=0.5)

    if plot_cbar:
        # Dynamically attach the colorbar to the drawn axes
        divider = make_axes_locatable(ax)
        orientation = "horizontal" if cbar_location in ["top", "bottom"] else "vertical"
        cax = divider.append_axes(cbar_location, size="3%", pad=0.5)

        sm = plt.cm.ScalarMappable(
            cmap=cmap, norm=plt.Normalize(vmin=obs_min, vmax=obs_max)
        )
        sm.set_array([])
        plt.colorbar(sm, cax=cax, orientation=orientation, label=cbar_label)

    format_ax(ax, extent=extent)
    return ax

# Convenience
def plot_side_by_side(
    field,
    title_field: str = "Velocity Field",
    title_obs: str = "Ray Paths",
    sources=sources,
    receivers=receivers,
    observations=data_obs,
    extent=model_grid_extent,
    figsize=(10, 8),
    field_cmap=cmc.batlow,
    obs_cmap=cmc.buda,
    cbar_location: Literal["left", "right", "top", "bottom"] = "bottom",
):
    """
    Generates a 1x2 subplot with shared axes for comparing the velocity field
    and the corresponding ray paths.
    """
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=figsize, constrained_layout=True, sharex=True, sharey=True
    )

    # Left Panel: Field
    plot_field(
        field, extent=extent, ax=ax1, cmap=field_cmap, cbar_location=cbar_location
    )
    plot_sender_receiver(sources=sources, receivers=receivers, extent=extent, ax=ax1)
    format_ax(ax1, title=title_field, extent=extent)
    add_legend(ax1)

    # Right Panel: Ray Paths
    plot_observations(
        sources=sources,
        receivers=receivers,
        observations=observations,
        extent=extent,
        ax=ax2,
        cmap=obs_cmap,
        cbar_location=cbar_location,
    )
    plot_sender_receiver(sources=sources, receivers=receivers, extent=extent, ax=ax2)
    format_ax(ax2, title=title_obs, extent=extent)

    return fig, (ax1, ax2)

# # Variogram and Covariance Model

# %%
# Estimate variogram - Cache to disk since it takes a few minutes to run
bin_centers_est, gamma_est = mem.cache(gs.vario_estimate)(
    (empirical_prior_model_grid_x, empirical_prior_model_grid_y),
    empirical_prior_samples,
    mesh_type="structured",
)

# %%
# Inspired by: https://geostat-framework.readthedocs.io/projects/gstools/en/v1.3.0/examples/03_variogram/01_find_best_model.html?highlight=variogram
models = {
    "Gaussian": gs.Gaussian,
    "Exponential": gs.Exponential,
    "Matern": gs.Matern,
    "Stable": gs.Stable,
    "Rational": gs.Rational,
    "Cubic": gs.Cubic,
    "Linear": gs.Linear,
    "Spherical": gs.Spherical,
    "Circular": gs.Circular,
    # "HyperSpherical": gs.HyperSpherical,
    # "SuperSpherical": gs.SuperSpherical,
    "JBessel": gs.JBessel,
}

fits = []
for model_name, model in models.items():
    fit_model = model(dim=2)
    # Note: gtol _very_ important as otherwise it does not go to convergence since residuals are just natively very small here
    para, pcov, r2 = fit_model.fit_variogram(
        bin_centers_est, gamma_est, return_r2=True, curve_fit_kwargs={"gtol": 1e-14}
    )
    fits.append(
        {
            "model_name": model_name,
            "model_class": model,
            "fit_model": fit_model,
            "para": para,
            "r2": r2,
        }
    )

fits = sorted(fits, key=lambda x: x["r2"], reverse=True)

# x_max = np.max(bin_centers_est)
x_max = 10.0
colors = cmc.batlowKS
styles = itertools.cycle(["--", "-", "-."])
for i, fit in enumerate(fits):
    style = next(styles)
    fit["i"] = i
    fit["style"] = style
    fit["color"] = colors(i + 3)

# Select model plot
select_fit = next(fit for fit in fits if fit["model_name"] == "Spherical")


# %%
covariance_model = select_fit["fit_model"]

# %% [markdown]
# ## Covariance Matrix

# %%
distance_matrix = get_distance_matrix(model_grid_x, model_grid_y)
cov_matrix = covariance_model.covariance(distance_matrix)

# Crucial: Add nugget to covariance matrix diagonal!
cov_matrix += np.identity(cov_matrix.shape[0]) * covariance_model.nugget

# %%
# Note: LLM
Nx = len(model_grid_x)
Ny = len(model_grid_y)

# 1. Find the flattened index of the center point
# Under "ij" indexing, y varies fastest, so index = i * Ny + j
center_i = Nx // 2
center_j = Ny // 2
center_idx = center_i * Ny + center_j

# 2. Extract the row and reshape back to (Nx, Ny)
dist_field = distance_matrix[center_idx, :].reshape(Nx, Ny)
cov_field = cov_matrix[center_idx, :].reshape(Nx, Ny)

# Physical extent for the plot axes
extent = [model_grid_x[0], model_grid_x[-1], model_grid_y[0], model_grid_y[-1]]

# %% [markdown]
# ## Prior Distribution

# %%
prior_dist = MultivariateNormalCholesky.from_covariance(
    cov=cov_matrix, mean=empirical_prior_mean, shape=model_grid_shape
)


# %%
# Note: No need to keep big arrays on GPU around, clean up
del distance_matrix
del cov_matrix

# %% [markdown]
# # Forward Model

# %%
# forward_model_solver = "skfmm:fmm"
# forward_model_solver_kwargs = {
#     # These arguments are tuned to work well for my system, depends on GPU, CPU, and memory bandwidth.
#     "chunk_size": 16,  # Trade-off between memory use and parallelization overhead
#     "parallel_args": {
#         "n_jobs": 4,
#     }
# }

forward_model_solver = "nanopinv:fsm"
forward_model_solver_kwargs = {
    "order": 2,
    "max_iter": 10_000,
    "tolerance": 1e-5,
}

forward_model = build_travel_time_points(
    sources,
    receivers,
    model_grid_x,
    model_grid_y,
    solver=forward_model_solver,
    solver_kwargs=forward_model_solver_kwargs,
)


# %%
# Check that forward model gives reasonable output for samples from the prior
forward_model_sanity_samples_model = prior_dist(skey(), num_samples=3)
forward_model_sanity_samples_data = forward_model(forward_model_sanity_samples_model)

# %% [markdown]
# # Likelihood
#

# %%
# Normally we would allow data_std to be an array as well, but here we know we have fixed STD for all observations,
# so we can make this constant and gain a bit of efficiency
@jax.jit(static_argnames=("data_std",))
def log_likelihood_gaussian(data, data_obs, data_std: float) -> Float[Array, ""]:
    """
    Compute the log-likelihood of observed data under a Gaussian noise model.
    """
    normalised_residual = (data_obs - data) / data_std
    log_likelihood = -0.5 * jnp.sum(normalised_residual**2)
    return jnp.sum(log_likelihood)


# %% [markdown]
# # Sample: EMC 1

# %%


# # %% [markdown]
# # # Sampler

# # %%
# # Configuration
# N_chains = 10
# keep_interval = 10

# init_step_size_min = 0.1
# init_step_size_max = 1.0

# betas_decay_base = 2.0

# N_steps_pre_burn_in = 500

# N_steps_tune_step_size_1 = 1500
# tune_interval_tune_step_size_1 = 300

# N_steps_tune_step_size_2 = 2000
# tune_interval_tune_step_size_2 = 400

# N_steps_tune_beta_1 = 3000
# tune_interval_beta_1 = 500

# target_chain_acceptance_rate = 0.25

# learning_rate_step_size = 1.0
# learning_rate_step_size_decay = 0.5
# learning_rate_beta = 1.0
# learning_rate_beta_decay = 0.9

# # %%
# # Temperature ladder for parallel tempering
# betas = initialize_betas(N_chains, base=betas_decay_base, last_is_zero=False)
# # Initial step size ladder
# step_sizes = np.linspace(init_step_size_min, init_step_size_max, N_chains)

# # %%
# proposal_dist = ProposalDistribution(dist=prior_dist, step_size=step_sizes)
# chain_template = ExtendedMetropolisChain(
#     beta=1.0,
#     proposal_dist=proposal_dist,
#     log_likelihood_fn=log_likelihood_gaussian,
#     forward_model=forward_model,
# )

# chain_spec = make_pytree_spec(
#     chain_template,
#     {
#         "beta": 0,  # Vectorize across chains
#         "proposal_dist.step_size": 0,  # Vectorize across chains
#         "*": None,  # Don't vectorize other fields
#     },
# )


# @eqx.filter_jit
# @eqx.filter_vmap(in_axes=(0, 0), out_axes=chain_spec)
# def make_chain(beta, step):
#     return ExtendedMetropolisChain(
#         beta=beta,
#         proposal_dist=ProposalDistribution(dist=prior_dist, step_size=step),
#         forward_model=forward_model,
#         log_likelihood_fn=log_likelihood_gaussian,
#     )


# chains = make_chain(betas, step_sizes)


# @eqx.filter_jit
# @eqx.filter_vmap(in_axes=(0, chain_spec, None))
# def init_state(k, c, obs_):
#     return c.get_iteration_state(prior_dist(k), obs_)

# # %%
# # Pre burn-in
# iter_states0 = init_state(skey(n=N_chains), chains, obs)

# pt = ParallelTemperingSampler(chains=chains, chain_axes_spec=chain_spec)

# # Pre-burn-in PT
# iter_states_burnt_in, pre_burn_in_history = pt.step_n(
#     N_steps_pre_burn_in,
#     key=skey(),
#     iter_states=iter_states0,
#     observations=obs,
#     keep_interval=keep_interval,
#     progress=True,
# )

# # %%
# pre_burn_in_history.plot_diagnostics()

# # %%
# # # Tune step-sizes only
# pt_tuned_1, iter_states_tuned_1, tune_history_1 = pt.tune_step_sizes(
#     n_steps_tune=N_steps_tune_step_size_1,
#     tune_interval=tune_interval_tune_step_size_1,
#     key=skey(),
#     iter_states=iter_states_burnt_in,
#     observations=obs,
#     target_chain_acceptance_rate=target_chain_acceptance_rate,
#     learning_rate=learning_rate_step_size,
#     learning_rate_decay=learning_rate_step_size_decay,
#     keep_interval=keep_interval,
#     progress=True,
# )


# # %%
# tune_history_1.plot_diagnostics()

# # %%
# # Tune temperatures
# # Try tune temperatures
# pt_tuned_betas_1, iter_states_tuned_betas_1, tune_history_betas_1 = (
#     pt_tuned_1.tune_betas(
#         n_steps_tune=N_steps_tune_beta_1,
#         tune_interval=tune_interval_beta_1,
#         key=skey(),
#         iter_states=iter_states_tuned_1,
#         observations=obs,
#         learning_rate=learning_rate_beta,
#         learning_rate_decay=learning_rate_beta_decay,
#         keep_interval=keep_interval,
#         method="czyz",
#         progress=True,
#     )
# )


# # %%
# tune_history_betas_1.plot_diagnostics(
#     window_chain=100, window_swap=tune_interval_beta_1
# )

# # %%
# # # Tune step-sizes only 2
# pt_tuned_step_sizes_2, iter_states_tuned_step_sizes_2, tune_history_step_sizes_2 = (
#     pt.tune_step_sizes(
#         n_steps_tune=N_steps_tune_step_size_2,
#         tune_interval=tune_interval_tune_step_size_2,
#         key=skey(),
#         iter_states=iter_states_tuned_betas_1,
#         observations=obs,
#         target_chain_acceptance_rate=target_chain_acceptance_rate,
#         learning_rate=learning_rate_step_size,
#         learning_rate_decay=learning_rate_step_size_decay,
#         keep_interval=keep_interval,
#         progress=True,
#     )
# )
# tune_history_step_sizes_2.plot_diagnostics(
#     window_chain=tune_interval_tune_step_size_2, window_swap=tune_interval_beta_1
# )

# # %%
# tune_history_step_sizes_2.plot_diagnostics(
#     window_chain=tune_interval_tune_step_size_2, window_swap=tune_interval_beta_1
# )

# # %%
# # Short unperturbed run to check autocorrelation
# iter_states_prod_1, prod_history_1 = pt_tuned_step_sizes_2.step_n(
#     n=2500,
#     key=skey(),
#     iter_states=iter_states_tuned_step_sizes_2,
#     observations=obs,
#     keep_interval=keep_interval,
#     progress=True,
# )
# tune_history_step_sizes_2.plot_diagnostics(
#     window_chain=tune_interval_tune_step_size_2, window_swap=tune_interval_beta_1
# )


# # %%
# keep_interval_prod_2 = 10
# # Production run
# iter_states_prod_2, prod_history_2 = pt_tuned_step_sizes_2.step_n(
#     n=10_000,
#     key=skey(),
#     iter_states=iter_states_prod_1,
#     observations=obs,
#     keep_interval=keep_interval_prod_2,
#     progress=True,
# )


# # %%
# prod_history_2.plot_diagnostics(
#     window_chain=keep_interval_prod_2, window_swap=tune_interval_beta_1, max_lag=200
# )

# # %%
# # Get accepted states from cold chain (beta=1.0)
# posterior_samples = prod_history_2.get_flat_cold_accepted_states(min_interval=25)
# posterior_samples.shape

# # %% [markdown]
# # # Posterior

# # %%
# print(posterior_samples.shape)

# plot_side_by_side(
#     posterior_samples.mean(axis=0),
#     title_field="Posterior Mean",
#     title_obs="Posterior Ray Paths",
#     field_cmap=cmc.batlow,
#     obs_cmap=cmc.buda,
# )

# ax = plot_field(posterior_samples.var(axis=0))
# format_ax(ax=ax, title="Posterior Variance")


# # %%
# # Check that forward of posterior matches observations well
# forward_posterior_samples = jax.vmap(forward_model)(posterior_samples)
# forward_posterior_samples_mean = forward_posterior_samples.mean(axis=0)
# forward_posterior_samples_p2_5 = jnp.percentile(forward_posterior_samples, 2.5, axis=0)
# forward_posterior_samples_p97_5 = jnp.percentile(
#     forward_posterior_samples, 97.5, axis=0
# )

# # And get some prior samples too
# prior_samples = prior_dist(skey(), num_samples=len(posterior_samples))
# forward_prior_samples = jax.vmap(forward_model)(prior_samples)
# forward_prior_samples_mean = forward_prior_samples.mean(axis=0)
# forward_prior_samples_p2_5 = jnp.percentile(forward_prior_samples, 2.5, axis=0)
# forward_prior_samples_p97_5 = jnp.percentile(forward_prior_samples, 95.5, axis=0)


# # %%
# plt.figure(figsize=(A4_WIDTH_IN, 5))
# i = jnp.arange(len(data_obs))

# cm = cmc.batlowKS
# alpha = 0.3

# # Plot prior
# plt.fill_between(
#     i,
#     forward_prior_samples_p2_5,
#     forward_prior_samples_p97_5,
#     color=cm(2),
#     alpha=0.3,
#     label="Prior 95% CI",
# )

# # Plot posterior
# plt.fill_between(
#     i,
#     forward_posterior_samples_p2_5,
#     forward_posterior_samples_p97_5,
#     color=cm(1),
#     alpha=0.3,
#     label="Posterior 95% CI",
# )

# plt.plot(i, forward_prior_samples_mean, "-", label="Prior Mean Prediction", color=cm(2))
# plt.plot(
#     i,
#     forward_posterior_samples_mean,
#     "-",
#     label="Posterior Mean Prediction",
#     color=cm(1),
# )

# data_obs_p2_5 = data_obs - 1.96 * obs.data_std
# data_obs_p97_5 = data_obs + 1.96 * obs.data_std
# plt.fill_between(
#     i,
#     data_obs_p2_5,
#     data_obs_p97_5,
#     color=cm(0),
#     alpha=0.3,
#     label="Observed Data 95% CI",
# )

# plt.plot(i, data_obs, "--", label="Observed Data", color=cm(0))
# plt.legend()
# plt.xlabel("Source-Receiver Pair Index, $i$")
# plt.ylabel("Travel Time [ns]")
# plt.title("Predictions vs Observed Data")
# plt.show()
# ...

# # %%
# # Residuals
# residuals_prior = forward_prior_samples - data_obs
# residuals_posterior = forward_posterior_samples - data_obs


# sns.histplot(
#     residuals_prior.flatten(),
#     color=cm(2),
#     label="Prior Residuals",
#     kde=True,
#     stat="density",
# )
# sns.histplot(
#     residuals_posterior.flatten(),
#     color=cm(1),
#     label="Posterior Residuals",
#     kde=True,
#     stat="density",
# )

# # Plot observation noise distribution
# x = np.linspace(-10, 10, 1000)
# obs_noise_dist = (1 / (obs.data_std * np.sqrt(2 * np.pi))) * np.exp(
#     -0.5 * (x / obs.data_std) ** 2
# )
# plt.plot(x, obs_noise_dist, color=cm(0), label="Observation Noise PDF")

# # Plot means
# plt.axvline(
#     residuals_prior.mean(), color=cm(2), linestyle="--", label="Prior Residual Mean"
# )
# plt.axvline(
#     residuals_posterior.mean(),
#     color=cm(1),
#     linestyle="--",
#     label="Posterior Residual Mean",
# )
# plt.axvline(0, color="k", linestyle=":", label="Zero Residual")

# plt.legend(frameon=True)
# plt.xlabel("Residual (Predicted - Observed) [ns]")
# plt.title("Residuals of Prior and Posterior Predictions")

# # %% [markdown]
# # # Porosity and Total Water In Place

# # %%
# # Water table definition
# elevations = model_grid_y
# water_table_elevation = 844.0
# water_table_y = model_grid_y[model_grid_y < water_table_elevation]
# extent_water_table = [
#     np.min(model_grid_x),
#     np.max(model_grid_x),
#     np.min(water_table_y),
#     np.max(water_table_y),
# ]


# @eqx.filter_jit
# def porosity(velocity):
#     # Equation (1) from assignment
#     return jnp.exp(-41.7 * velocity + 2.03)


# def extract_water_table(field, elevations, cut_off_elevation: float):
#     mask = elevations < cut_off_elevation
#     water_table = field[:, mask]
#     return water_table


# mean_velocity_field = posterior_samples.mean(axis=0)
# mean_porosity_field = porosity(mean_velocity_field)

# # Plot porosity
# plot_side_by_side(
#     mean_porosity_field,
#     title_field="Mean Porosity Field",
# )

# # Plot porosity in water table
# mean_porosity_water_table = extract_water_table(
#     mean_porosity_field, elevations, water_table_elevation
# )

# plot_side_by_side(
#     mean_porosity_water_table,
#     title_field="Mean Porosity in Water Table",
#     title_obs="Ray Paths in Water Table",
#     extent=extent_water_table,
# )


# # Calculate TWIP
# def total_water_in_place(porosity_field, model_grid_dx, model_grid_dy):
#     total_porosity = jnp.sum(porosity_field)
#     cell_area = model_grid_dx * model_grid_dy
#     water_volume = total_porosity * cell_area
#     return water_volume


# total_water = total_water_in_place(
#     mean_porosity_water_table, model_grid_dx, model_grid_dy
# )
# print(f"Total water in place: {total_water:.2f} m^3")

# # %%
# # Try to do covmodel on posterior, just to check

# posterior_mean = posterior_samples.mean()
# print(f"{posterior_mean}")

# covmodel_posterior = gs.Spherical(dim=2)
# bin_centers_est_posterior, gamma_est_posterior = mem.cache(gs.vario_estimate)(
#     (model_grid_x, model_grid_y),
#     posterior_samples[:250],
#     mesh_type="structured",
# )

# covmodel_posterior.fit_variogram(
#     bin_centers_est_posterior, gamma_est_posterior, curve_fit_kwargs={"gtol": 1e-14}
# )

# # Plot fits
# plt.figure()
# plt.scatter(
#     bin_centers_est_posterior,
#     gamma_est_posterior,
#     color="k",
#     label="Empirical Posterior",
# )
# covmodel_posterior.plot(
#     ax=plt.gca(), x_max=x_max, label="Fitted Posterior Covariance Model"
# )
# plt.xlabel("Distance, $h$ [m]")
# plt.ylabel("Semivariance, $γ(h)$ [ns$^2$]")
# plt.title("Empirical variogram of posterior samples and fitted model")

# print(f"{covmodel_posterior}")

# # covmodel_posterior.length_scale = 2

# # Now try and generate one to have a test sample to use
# distance_matrix_posterior = get_distance_matrix(model_grid_x, model_grid_y)
# cov_matrix_posterior = covmodel_posterior.covariance(distance_matrix_posterior)
# cov_matrix_posterior += (
#     np.identity(cov_matrix_posterior.shape[0]) * covmodel_posterior.nugget
# )
# dist_obj_posterior = MultivariateNormalCholesky.from_covariance(
#     cov=cov_matrix_posterior,
#     mean=posterior_samples.mean(axis=0),
#     shape=model_grid_shape,
# )

# dist_posterior_samples = dist_obj_posterior(skey(), num_samples=3)
# # Plot 3 samples
# for i in range(3):
#     plot_field(dist_posterior_samples[i])

# # %%
# # Calculate TWIP across all
# water_table_mask = model_grid_y < water_table_elevation

# accepted_velocity_fields_water_table = posterior_samples[:, :, water_table_mask]
# print(accepted_velocity_fields_water_table.shape)


# @eqx.filter_vmap(in_axes=(0, None, None))
# @jax.jit(
#     static_argnames=(
#         "model_grid_dx",
#         "model_grid_dy",
#     )
# )
# def total_water_in_place_from_velocity_field(
#     velocity_field, model_grid_dx, model_grid_dy
# ) -> float:
#     porosity_water_table = porosity(velocity_field)
#     twip = total_water_in_place(porosity_water_table, model_grid_dx, model_grid_dy)

#     return twip


# twip_posterior = total_water_in_place_from_velocity_field(
#     accepted_velocity_fields_water_table, model_grid_dx, model_grid_dy
# )
# display(twip_posterior)
# twip_mean = twip_posterior.mean()
# twip_sd = twip_posterior.std()
# print(f"TWIP_posterior = {twip_mean:.2f} ± {twip_sd:.2f}")

# display(
#     total_water_in_place_from_velocity_field(
#         posterior_samples, model_grid_dx, model_grid_dy
#     )
# )

# sns.histplot(twip_posterior, edgecolor="k")

# # %%
# # Test the prior distribution
# prior_velocity_fields = prior_dist(skey(), num_samples=100)
# prior_velocity_fields_water_table = prior_velocity_fields[:, :, water_table_mask]
# twip_prior = total_water_in_place_from_velocity_field(
#     prior_velocity_fields_water_table, model_grid_dx, model_grid_dy
# )

# display(twip_prior)
# twip_mean = twip_prior.mean()
# twip_sd = twip_prior.std()

# print(f"TWIP_prior = {twip_mean:.2f} ± {twip_sd:.2f}")
# display(
#     total_water_in_place_from_velocity_field(
#         prior_velocity_fields, model_grid_dx, model_grid_dy
#     )
# )

# # %%
# sns.histplot(twip_prior, edgecolor="k", label="Prior")

# # %%
# accepted_states.shape


