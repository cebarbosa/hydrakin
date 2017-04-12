# -*- coding: utf-8 -*-
"""

Created on 10/04/2017

@Author: Carlos Eduardo Barbosa

Zoomed maps for HCC 007

"""
import os

import numpy as np
import matplotlib.pyplot as plt

from config import *
import canvas as cv

def get_positions_by_slits(slits, canvas):
    """ Matches two different tables using the spectra column. """
    xy = []
    for i, slit in enumerate(slits):
        index = canvas.slits.ids.index(slit)
        xy.append([canvas.slits.x[index], canvas.slits.y[index]])
    return np.array(xy)

def find_chart_hcc007(intable):
    """ Produces maps for HCC 007."""
    specs = np.loadtxt(intable, usecols=(0,), dtype=str)
    ids = [x.split("n3311")[1].split(".")[0] for x in specs]
    canvas = cv.CanvasImage("vband")
    fig, ax = plt.subplots(1, 1, figsize=(3.54, 2.5))
    ax.minorticks_on()
    ax.set_aspect("equal")
    xy = get_positions_by_slits(ids, canvas)
    canvas.make_contours(lw=0.3)
    canvas.draw_slits(ax, slit_type=3)
    ax.set_xlim(-2.5,-10.5)
    ax.set_ylim(-31.5,-26)
    for a, b, c in zip(xy[:,0], xy[:,1], ids):
        ax.text(a, b, c, weight="bold", color="C0", fontsize=8)
    ax.set_xlabel("X (kpc)")
    ax.set_ylabel("Y (kpc)")
    plt.subplots_adjust(right=0.98, top=0.99, bottom=0.135)

    plt.savefig(os.path.join(figures_dir, "find_chart_hcc007.png"),
                dpi=250)
    return

if __name__ == "__main__":
    plt.style.use("seaborn-paper")
    table = os.path.join(results_dir, "ppxf_results_hcc007.dat")
    find_chart_hcc007(table)