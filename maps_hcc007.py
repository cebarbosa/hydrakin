# -*- coding: utf-8 -*-
"""

Created on 10/04/2017

@Author: Carlos Eduardo Barbosa

Zoomed maps for HCC 007

"""
import os
import fileinput

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize

from config import *
import canvas as cv
from newcolorbars import cubelaw

def get_positions_by_slits(slits, canvas):
    """ Matches two different tables using the spectra column. """
    xy = []
    for i, slit in enumerate(slits):
        index = canvas.slits.ids.index(slit)
        xy.append([canvas.slits.x[index], canvas.slits.y[index]])
    return np.array(xy)

def find_chart_hcc007(intable):
    """ Produces finding chart for HCC 007."""
    specs = np.loadtxt(intable, usecols=(0,), dtype=str)
    ids = [x.split("n3311")[1].split(".")[0] for x in specs]
    canvas = cv.CanvasImage("vband")
    fig, ax = plt.subplots(1, 1, figsize=(3.54, 2.5))
    ax.minorticks_on()
    ax.set_aspect("equal")
    xy = get_positions_by_slits(ids, canvas)
    canvas.make_contours(lw=0.3)
    canvas.draw_slits(ax, slit_type=3)
    canvas.draw_slits(ax, slit_type=1)
    ax.set_xlim(-2.5,-10.5)
    ax.set_ylim(-31.5,-26)
    for a, b, c in zip(xy[:,0], xy[:,1], ids):
        ax.text(a, b, c, weight="bold", color="C0", fontsize=8)
    ax.set_xlabel("X (kpc)")
    ax.set_ylabel("Y (kpc)")
    plt.subplots_adjust(right=0.98, top=0.99, bottom=0.135)
    plt.savefig(os.path.join(figures_dir, "find_chart_hcc007.png"),
                dpi=250)
    plt.show()
    return

def make_kinematics():
    """ Produces maps for the kinematics. """
    tables =  [os.path.join(results_dir, "ppxf_results_hcc007.dat"),
               os.path.join(results_dir, "ppxf_results_best.dat")]
    specs = np.loadtxt(fileinput.input(tables), usecols=(0,), dtype=str)
    data = np.loadtxt(fileinput.input(tables), usecols=(1,3,5,7)).T
    ids = [x.split("n3311")[1].split(".")[0] for x in specs]
    canvas = cv.CanvasImage("vband")
    idx = [canvas.slits.ids.index(iid) for iid in ids]
    rects = canvas.calc_vertices(canvas.slits.x[idx], canvas.slits.y[idx],
                               canvas.slits.w[idx], canvas.slits.l[idx],
                               canvas.slits.ang[idx] + canvas.posangle)
    xy = get_positions_by_slits(ids, canvas)
    cmap = "Spectral_r"
    names = ["vel", "sigma", "h3", "h4"]
    cb_fmts = ["%i", "%i", "%.2f", "%.2f"]
    zlims = [[4500, 4870], [50, 110], [-0.1, .1], [-0.05, .1]]
    labels = ["V (km/s)", "$\sigma$ (km/s)", "$h_3$", "$h_4$"]
    for i, d in enumerate(data):
        fig, ax = plt.subplots(1, 1, figsize=(3.54, 3))
        ax.minorticks_on()
        ax.set_aspect("equal")
        ax.imshow(canvas.data, extent = canvas.extent, cmap="bone",
                  vmin=18, vmax=24, origin="bottom")
        norm = Normalize(vmin=zlims[i][0], vmax=zlims[i][1])
        canvas.make_contours(lw=0.3)
        coll = PolyCollection(rects, array=d,
                              cmap=cmap, edgecolors='w', norm=norm,
                              linewidths=0.4)
        ax.add_collection(coll)
        ax.set_xlim(0., -12.)
        ax.set_ylim(-34, -24)
        ax.set_xlabel("X (kpc)")
        ax.set_ylabel("Y (kpc)")
        cbar_pos = [0.65, 0.25, 0.25, 0.05]
        cbaxes = fig.add_axes(cbar_pos)
        cbar = plt.colorbar(coll, cax=cbaxes, orientation="horizontal",
                            format=cb_fmts[i])
        cbar.set_ticks(np.linspace(zlims[i][0], zlims[i][1], 4))
        cbar.ax.xaxis.set_label_position('bottom')
        cbar.ax.xaxis.set_ticks_position('bottom')
        cbar.set_label(labels[i], labelpad = -34)
        cl = plt.getp(cbar.ax, 'ymajorticklabels')
        ax = cbar.ax
        plt.subplots_adjust(right=0.97, top=0.97, bottom=0.135)
        plt.savefig(os.path.join(figures_dir,
                    "hcc007_{}.png".format(names[i])), dpi=250)
    return


if __name__ == "__main__":
    plt.style.use("seaborn-paper")
    table = os.path.join(results_dir, "ppxf_results_hcc007.dat")
    find_chart_hcc007(table)
    # make_kinematics()