import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import chaosmagpy as cp


RADIUS_EARTH = 6371.2  # Earth radius in km
RADIUS_CORE = 3480.0  # Core radius in km


# Given in assignment
# Global plot with polar areas
def plot_global(
    lon, lat, data, point_size, title="", cbar_label="", cmap=None, limit=None
):
    if limit == None:
        limit = np.max(abs(data))
    else:
        limit = np.max(abs(np.array(limit)))
    # create figure
    fig = plt.figure(figsize=(9, 5))
    # make array of axes
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
    # Iterate over axes
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
            linewidth=0.5,
            ylocs=np.linspace(-90, 90, num=7),  # parallels
            xlocs=np.linspace(-180, 180, num=13),
            color="grey",
            alpha=0.6,
            linestyle="-",
        )  # meridians
        ax.coastlines(linewidth=0.5)
    # Add colorbar
    # inset axes into global map and move upwards
    cax = inset_axes(
        axes[-1], width="55%", height="10%", loc="upper center", borderpad=-9
    )
    # use last artist for the colorbar
    clb = plt.colorbar(pc, cax=cax, extend="both", orientation="horizontal")
    clb.set_label("{}".format(cbar_label), fontsize=12)

    # Title
    plt.suptitle("{}".format(title))
    # Adjust plot
    plt.subplots_adjust(
        top=0.985, bottom=0.015, left=0.008, right=0.992, hspace=0.0, wspace=0.0
    )
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