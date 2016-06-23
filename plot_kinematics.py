# -*- coding: utf-8 -*-
"""

Created on 21/06/16

@author: Carlos Eduardo Barbosa

"""
import os

import numpy as np
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

from config import *
import canvas as cv

def get_ventimiglia2010():
    """ Read data from Ventimiglia+ 2010. """
    ra, dec, v, verr, sig, sigerr= np.loadtxt(os.path.join(tables_dir,
                                "ventimiglia_kin.txt"), usecols=np.arange(6)).T
    x, y = radec2xy(ra, dec)
    r, theta = xy2polar(x, y)
    velocity_offset = 72.0
    z = np.nan * np.zeros_like(v)
    return np.column_stack((x, y, r, theta, v+velocity_offset,verr, sig,
                            sigerr, z, z, z, z))

def get_richtler():
    """ Retrieve data from Richtler et al. 2011. """
    files = [os.path.join(tables_dir, "results_n3311maj_sn30.txt"),
             os.path.join(tables_dir, "results_n3311min_sn30.txt")]
    tables = []
    for i, fname in enumerate(files):
        ras, v, verr, sig, sigerr = np.loadtxt(fname, usecols=(1,6,8,7,9)).T
        r = np.sign(ras) * D * 1000 * np.tan(np.deg2rad(np.abs(ras)/3600))
        if i == 0:
            theta = np.where(r> 0, 40 * np.ones_like(r),
                           (40) * np.ones_like(r))
        else:
            theta = np.where(r > 0, 115 * np.ones_like(r),
                            (115.) * np.ones_like(r))
        x = r * np.sin(np.deg2rad(theta))
        y = r * np.cos(np.deg2rad(theta))
        r, theta = xy2polar(x, y)
        velocity_offset = -14.0
        z = np.nan * np.zeros_like(v)
        table = np.column_stack((x, y, r, theta, v + velocity_offset, verr,
                                 sig, sigerr, z, z, z, z))
        if i == 0:
            tables = table
        else:
            tables = np.vstack((tables, table))
    return tables

def radec2xy(ra, dec):
    """ Convert from coordinates to X and Y relative to the center of N3311."""
    x = D * 1000 * np.tan(np.deg2rad(ra - ra0))
    y = D * 1000 * np.tan(np.deg2rad(dec - dec0))
    return x, y

def xy2polar(x, y):
    """ Convert from cartesian to polar coordinates """
    r = np.sqrt(x**2 + y**2)
    theta = np.rad2deg(np.arctan2(x,y))
    theta[theta < 0] += 360.
    return r, theta


def plot_cones(datasets, pas=None, dpa=22.5):
    """ Make radial plots using conic sections in the XY plane. """
    # General setting for the plots
    colors = ("r", "b", "g")
    symbols = ("o", "s", "^")
    ylims = [[3300, 4600], [100, 600], [-.3, .3], [-.3,.3]]
    xlims = [-40, 40]
    ylabels = [r"$V_{\rm{LOS}}$ (km/s)", r"$\sigma_{\rm{LOS}}$ (km/s)",
               r"$h_3$", r"$h_4$"]
    mec = "0.7" # Grayscale level for bars and symbol edge colors
    names = ["rad_vel", "rad_sigma", "rad_h3", "rad_h4"]
    sn_min = [15, 15, 30, 30]
    sn = np.loadtxt(intable, usecols=(14,))
    ##########################################################################
    # Set the default position angles
    if pas is None:
        pas = np.linspace(0, 180, 5)[:-1] + pa0
    ##########################################################################
    for mm in range(4):
        for i, d in enumerate(datasets):
            x, y, r, theta = d[:,:4].T
            moment = d[:,np.arange(4,12,2)[mm]].T
            error = np.clip(d[:,np.arange(5,13,2)[mm]].T, 0, 1000)
            fig = plt.figure(1, figsize=(6, 8.5))
            for j, pa in enumerate(pas):
                ###############################################################
                # Select conic regions
                idx1, idx2 = [], []
                for a in np.linspace(-360,360,3):
                    idx1 += np.argwhere((theta > pa - dpa + a) &
                                        (theta < pa + dpa + a)).tolist()
                    idx2 += np.argwhere((theta > pa - dpa + a + 180) &
                                        (theta < pa + dpa + a + 180)).tolist()
                    #==========================================================
                    # Include central points of our dataset in all subplots
                    if i == 0:
                        idx1 += np.argwhere((theta > pa - 90 + a) &
                                            (theta < pa + 90 + a) &
                                            (r < 8)).tolist()
                        idx2 += np.argwhere((theta > pa - 90 + a + 180) &
                                            (theta < pa + 90 + a + 180) &
                                            (r < 8)).tolist()
                    #==========================================================
                idx1 = np.unique(np.array(idx1).ravel())
                idx2 = np.unique(np.array(idx2).ravel())
                #=============================================================
                # S/N cut for our dataset
                if i == 0:
                    idxsn = np.where(sn > sn_min[mm])
                    idx1 = np.intersect1d(idx1, idxsn)
                    idx2 = np.intersect1d(idx2, idxsn)
                #=============================================================
                ##############################################################
                # Produces figure
                ax = plt.subplot(len(pas), 1, j+1)
                ax.minorticks_on()
                if len(idx1) > 0:
                    ax.errorbar(r[idx1], moment[idx1], yerr=error[idx1],
                        fmt=symbols[i], ecolor=mec, c=colors[i],
                        mec=mec, ms=8, zorder=-i)
                if len(idx2) > 0:
                    ax.errorbar(-r[idx2], moment[idx2], yerr=error[idx2],
                        fmt=symbols[i], ecolor=mec, c=colors[i],
                        mec=mec, ms=8, zorder=-i)
                ax.set_xlim(xlims)
                ax.set_ylim(ylims[mm])
                ax.set_ylabel(ylabels[mm])
                ax.axvline(x=0, ls="--", c="k")
                ax.annotate("PA={0:.1f}$\pm${1:.1f}".format(pa, dpa),
                            xy=(0.5, 0.8), xycoords='axes fraction',
                            fontsize=14, horizontalalignment='center',
                            verticalalignment='bottom',
                            bbox=dict(boxstyle="round, pad=0.3", fc="w"))
                if mm > 1:
                    ax.axhline(y=0, ls="--", c="k")
            ax.set_xlabel("R (kpc)")
        plt.subplots_adjust(left=0.14, right=0.96, bottom=0.07, top=0.98)
        plt.savefig(os.path.join(figures_dir, "{0}.png".format(names[mm])))
        plt.clf()
    return

def plot_rings(datasets, radius=None):
    """ Make azimuthal plots in radial sections. """
    # General setting for the plots
    colors = ("r", "b", "g")
    symbols = ("o", "s", "^")
    ylims = [[3300, 4600], [100, 700], [-.3, .3], [-.3,.3]]
    xlims = [0, 360]
    ylabels = [r"$V_{\rm{LOS}}$ (km/s)", r"$\sigma_{\rm{LOS}}$ (km/s)",
               r"$h_3$", r"$h_4$"]
    mec = "0.7" # Grayscale level for bars and symbol edge colors
    names = ["pa_vel", "pa_sigma", "pa_h3", "pa_h4"]
    sn_min = [10, 10, 10, 10]
    sn = np.loadtxt(intable, usecols=(14,))
    ##########################################################################
    # Set the default position angles
    if radius is None:
        radius = np.linspace(0,40,5)
    ##########################################################################
    for mm in range(4):
        for i, d in enumerate(datasets):
            x, y, r, theta = d[:,:4].T
            moment = d[:,np.arange(4,12,2)[mm]].T
            error = np.clip(d[:,np.arange(5,13,2)[mm]].T, 0, 1000)
            fig = plt.figure(1, figsize=(6, 8.5))
            for j in range(len(radius) - 1):
                rmin = radius[j]
                rmax = radius[j+1]
                ###############################################################
                # Select regions
                idx = np.argwhere((r >= rmin) & (r < rmax))
                # S/N cut for our dataset
                if i == 0:
                    idxsn = np.where(sn > sn_min[mm])
                    idx = np.intersect1d(idx, idxsn)
                ##############################################################
                idx = idx.ravel()
                # Produces figure
                ax = plt.subplot(len(radius)-1, 1, len(radius)-1-j)
                ax.minorticks_on()
                if len(idx) > 0:
                    ax.errorbar(theta[idx], moment[idx], yerr=error[idx],
                        fmt=symbols[i], ecolor=mec, c=colors[i],
                        mec=mec, ms=8, zorder=-i)
                ax.set_xlim(xlims)
                ax.set_ylim(ylims[mm])
                ax.set_ylabel(ylabels[mm])
                ax.axvline(x=0, ls="--", c="k")
                ax.annotate("$R\in [{0:.1f},{1:.1f}[$".format(rmin, rmax),
                            xy=(0.75, 0.8), xycoords='axes fraction',
                            fontsize=14, horizontalalignment='center',
                            verticalalignment='bottom',
                            bbox=dict(boxstyle="round, pad=0.3", fc="w"))
                if mm > 1:
                    ax.axhline(y=0, ls="--", c="k")
                if j == 0:
                    ax.set_xlabel("PA (degree)")
                ax.axvline(x=63, ls="--", c="k")
                ax.axvline(x=63+180, ls="--", c="k")
        plt.subplots_adjust(left=0.14, right=0.96, bottom=0.07, top=0.98)
        plt.savefig(os.path.join(figures_dir, "{0}.png".format(names[mm])))
        plt.clf()
    return

def cones_vertical(pas=None, dpa=22.5):
    """ Make map indicating the conic sections. """
    if pas is None:
        pas = np.linspace(0, 180, 5)[:-1] + pa0
    canvas = cv.CanvasImage("vband")
    canvas.data_smooth = ndimage.gaussian_filter(canvas.data, 3, order=0.)
    contours = np.linspace(19, 23.5, 10)
    extent = canvas.extent
    fig = plt.figure(1, figsize=(2.8, 8.5))
    l = 1000
    for j, pa in enumerate(pas):
        ax = plt.subplot(len(pas), 1, j+1)
        ax.minorticks_on()
        ax.contour(canvas.data_smooth, contours,
                   extent=extent, colors="k", linewidths=0.5)
        canvas.draw_slits(ax, fc="r", ec="r")
        canvas.draw_literature(ax)
        ax.set_xlim(40, -40)
        ax.set_ylim(-40, 40)
        ax.set_ylabel("Y (kpc)")
        ax.tick_params(labelsize=10)
        x1 = np.array([l * np.sin(np.deg2rad(pa + dpa)),
                       -l * np.sin(np.deg2rad(pa + dpa))])
        y1 = np.array( [l * np.cos(np.deg2rad(pa + dpa)),
                        - l * np.cos(np.deg2rad(pa + dpa))])
        x2 = np.array([l * np.sin(np.deg2rad(pa - dpa)),
                       - l * np.sin(np.deg2rad(pa - dpa))])
        y2 = np.array([l * np.cos(np.deg2rad(pa - dpa)),
                       - l * np.cos(np.deg2rad(pa - dpa))])
        line1 = interp1d(x1, y1)
        line2 = interp1d(x2, y2)
        x = np.linspace(-40, 40, 1000)
        if pa < 180.:
            ax.fill_between(x, line1(x), line2(x), color="0.5", alpha=0.6)
        else:
            ax.fill_between(x, line1(x), -line2(x), color="0.5", alpha=0.6)
    ax.set_xlabel("X (kpc)")
    plt.subplots_adjust(left=0.23, right=0.93, bottom=0.07, top=0.98)
    plt.savefig(os.path.join(figures_dir, "cones.png"))

if __name__ == "__main__":
    os.chdir(results_dir)
    np.set_printoptions(5, suppress=True)
    ###########################################################################
    # Loading data
    intable = "results.tab"
    data = np.loadtxt(intable, usecols=(1,2,3,4,5,6,7,8,9,10,11,12))
    idx = data[:,3] < 0
    data[idx,3] += 360.
    v10 = get_ventimiglia2010()
    r11 = get_richtler()
    ###########################################################################
    # Radial plots in conic sections
    # plot_cones((data, v10, r11), pas=None, dpa=22.5)
    ##########################################################################
    # Azimuthal plots
    plot_rings((data, v10, r11))
    ##########################################################################
    # cones_vertical()
