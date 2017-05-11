# -*- coding: utf-8 -*-
"""

Created on 11/05/2017

@Author: Carlos Eduardo Barbosa

Plots for the kinematics of HC 007.

"""
import os

import numpy as np
import matplotlib.pyplot as plt

from config import *
import canvas as cv


def plot_sma_losvd():
    """ Produces plot for LOSVD as a function of the radius. """
    # Reading table
    filename = os.path.join(results_dir, "ppxf_results_hcc007.dat")
    specs = np.loadtxt(filename, usecols=(0,), dtype=str)
    data = np.loadtxt(filename, usecols=(1,3,5,7))
    errors = np.loadtxt(filename, usecols=(2,4,6,8))
    ##########################################################################
    # Getting (x,y)
    canvas = cv.CanvasImage("vband")
    ids = [x.split("n3311")[1].replace(".fits", "") for x in specs]
    xy = []
    for slit in ids:
        idx = canvas.slits.ids.index(slit)
        xy.append([canvas.slits.x[idx], canvas.slits.y[idx]])
    x, y = np.array(xy).T
    ##########################################################################
    # Calculating isophotal distances
    idx0 = canvas.slits.ids.index("inn2_s39a")
    x0 = canvas.slits.x[idx0]
    y0 = canvas.slits.y[idx0]
    pa0 = -51
    q = 0.7
    sma = calc_isophotes(x, y, x0, y0, pa0, q)
    pa = np.rad2deg(np.arctan2(x - x0, y - y0))
    r = np.sqrt((x-x0)**2 + (y-y0)**2)
    print x, y
    sign = np.where(np.cos(np.radians(pa - pa0)) > 0, 1, -1)
    # Plotting
    plt.style.use("seaborn-paper")
    fig = plt.figure(figsize=(3.54, 5))
    ylabels = ["V (km s$^{\\rm -1}$)", "$\sigma$ (km s$^{\\rm -1}$)",
               "$h_3$", "$h_4$"]
    for i in range(4):
        ax = plt.subplot(4, 1,i+1)
        ax.minorticks_on()
        ax.errorbar(sign * sma,  data[:,i], yerr=errors[:,i], fmt="o",
                    ecolor="0.8")
        if i < 3:
            ax.xaxis.set_ticklabels([])
        else:
            ax.set_xlabel("semi-major axis (kpc)")
        ax.set_ylabel(ylabels[i])
    plt.subplots_adjust(left=0.15, bottom=0.09, top=0.98, right=0.98)
    print figures_dir
    plt.savefig(os.path.join(figures_dir, "sma_losvd.png"), dpi=250)

def calc_isophotes(x, y, x0, y0, PA, q):
    """ Calculate isophotes """
    x = np.copy(x) - x0
    y = np.copy(y) - y0
    shape = x.shape
    theta = np.radians(PA)
    c, s = np.cos(theta), np.sin(theta)
    rot = np.array([[s, -c], [c, s]])
    xy = np.dot(np.column_stack((x.flatten(), y.flatten())), rot).T
    x = np.reshape(xy[0], newshape=shape)
    y = np.reshape(xy[1], newshape=shape)
    return np.sqrt(np.power(x, 2) + np.power(y / q, 2))


if __name__ == "__main__":
    plot_sma_losvd()
