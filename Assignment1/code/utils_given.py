import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import chaosmagpy as cp


RADIUS_EARTH = 6371.2  # Earth radius in km
RADIUS_CORE = 3480.0  # Core radius in km

def plot_global(
    lon,
    lat,
    data,
    point_size,
    title="",
    cbar_label="",
    cmap=None,
    limit=None,
    half_page: bool = False,
):
    if limit is None:
        limit = np.max(np.abs(data))
    else:
        limit = np.max(np.abs(np.array(limit)))

    if half_page:
        fig = plt.figure(figsize=(4.0, 4.0 * 0.62))
        grid_lw = 0.35
        coast_lw = 0.35
        cax_width = "50%"
        cax_height = "12%"
        cax_borderpad = -5
        cbar_label_size = 8
        cbar_tick_size = 7
        title_size = 11
        adjust_kwargs = dict(
            top=0.90, bottom=0.05, left=0.03, right=0.97, hspace=0.02, wspace=0.02
        )
    else:
        fig = plt.figure(figsize=(9, 5))
        grid_lw = 0.5
        coast_lw = 0.5
        cax_width = "55%"
        cax_height = "10%"
        cax_borderpad = -9
        cbar_label_size = 12
        cbar_tick_size = None
        title_size = None
        adjust_kwargs = dict(
            top=0.985, bottom=0.015, left=0.008, right=0.992, hspace=0.0, wspace=0.0
        )

    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1], height_ratios=[0.35, 0.65])
    axes = []
    axes.append(
        plt.subplot(
            gs[0, 0], projection=ccrs.NearsidePerspective(central_latitude=90.0)
        )
    )
    axes.append(
        plt.subplot(
            gs[0, 2], projection=ccrs.NearsidePerspective(central_latitude=-90.0)
        )
    )
    axes.append(plt.subplot(gs[1, :], projection=ccrs.Mollweide()))

    for ax in axes:
        pc = ax.scatter(
            lon,
            lat,
            c=data,
            s=point_size,
            cmap=cmap,
            vmin=-limit,
            vmax=limit,
            transform=ccrs.PlateCarree(),
        )
        ax.gridlines(
            linewidth=grid_lw,
            ylocs=np.linspace(-90, 90, num=7),
            xlocs=np.linspace(-180, 180, num=13),
            color="grey",
            alpha=0.6,
            linestyle="-",
        )
        ax.coastlines(linewidth=coast_lw)

    cax = inset_axes(
        axes[-1],
        width=cax_width,
        height=cax_height,
        loc="upper center",
        borderpad=cax_borderpad,
    )
    clb = plt.colorbar(pc, cax=cax, extend="both", orientation="horizontal")
    clb.set_label("{}".format(cbar_label), fontsize=cbar_label_size)
    if cbar_tick_size is not None:
        clb.ax.tick_params(labelsize=cbar_tick_size)

    if title_size is None:
        plt.suptitle("{}".format(title))
    else:
        plt.suptitle("{}".format(title), fontsize=title_size)

    plt.subplots_adjust(**adjust_kwargs)
    return fig


def plot_power_spectrum(m, N: int, radius: float = RADIUS_CORE):
    # SH power spectrum and Br on grid at CMB
    ps = cp.model_utils.power_spectrum(m, radius=radius)

    # Plot Power spectrum (Unphysical divergence at high SH degree)
    n = np.arange(1, N + 1)
    plt.figure(figsize=(9, 4))
    plt.semilogy(n, ps)
    plt.xlabel("degree n")
    plt.ylabel(r"Power [nT$^2$]")
    plt.title(r"Power spectra of estimated geomagnetic field at CMB $r={}$km".format(radius))
    nmax = np.max([len(n)])
    plt.xticks(np.arange(1, nmax + 1))
    plt.grid()
    plt.show()