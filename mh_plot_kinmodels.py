# -*- coding: utf-8 -*-
"""

Created on 27/04/17

@author: Michael Hilker

"""
import os

import numpy as np
import matplotlib.pyplot as plt

from config import *
# from operator import itemgetter
import canvas as cv

def plot_models1():
    """ Produces plot for velocity dispersion as a function of radius. """
    # Selection of PA ranges
    radiusne, sigmane, esigne = seltabne[:,[0,1,2]].T
    radiussw, sigmasw, esigsw = seltabsw[:,[0,1,2]].T
    # Reading model table
    modelname = os.path.join(results_dir, "model_tmp.dat")
    radm, mod1 = np.loadtxt(modelname, usecols=(0, 1), unpack=True)
    # Plotting
    plt.style.use("seaborn-paper")
    fs = _small_fig_settings()
    fig = plt.figure(1, figsize=(fs["width"], fs["height"]))
    ax = plt.subplot(1, 1, 1)
    ax.minorticks_on()
    ax.set_xlim(0, 39)
    ax.set_ylim(0,950)
    ax.axvline(x=re, ls="--", c="0.5")
    ax.errorbar(radius, sigma, yerr=esig, fmt="ko", ecolor="0.8")
    # ax.errorbar(radiussw, sigmasw, yerr=esigsw, fmt="bo", ecolor="0.8")
    # ax.errorbar(radiusne, sigmane, yerr=esigne, fmt="ro", ecolor="0.8")
    # ax.plot(radm, mod1, 'go-')
    ax.set_xlabel(r"$R$ (kpc)", size=12)
    ax.set_ylabel(r"$\sigma_{\rm{LOS}}$ (km/s)", size=12)
    plt.subplots_adjust(left=0.14, right=fs["right"], bottom=0.18,
                        top=fs["top"], hspace=fs["hspace"])
    print figures_dir
    plt.savefig(os.path.join(figures_dir, "disp_profile_models.png"), dpi=250)

def plot_velsigma():
    """ Produces plot for V/sigma as a function of radius. """
    # Selection of PA ranges
    radiusne, vsigne, evsigne = seltabne[:,[0,3,4]].T
    radiussw, vsigsw, evsigsw = seltabsw[:,[0,3,4]].T
    # Fitting
    # fitall = np.polyfit(radius, vsig, 2, w=(1/evsig))
    fitne = np.polyfit(radiusne, vsigne, 2, w=(1/np.sqrt(evsigne)))
    fitsw = np.polyfit(radiussw, vsigsw, 2, w=(1/np.sqrt(evsigsw)))
    # Plotting
    plt.style.use("seaborn-paper")
    fs = _large_fig_settings()
    fig = plt.figure(1, figsize=(fs["width"], fs["height"]))
    ax = plt.subplot(1, 1, 1)
    ax.minorticks_on()
    ax.set_xlim(0, 39)
    ax.set_ylim(-0.2, 2.4)
    ax.axvline(x=re, ls="--", c="0.5")
    ax.errorbar(radius, vsig, yerr=evsig, fmt="ko", ecolor="0.8")
    ax.errorbar(radiussw, vsigsw, yerr=evsigsw, fmt="bo", ecolor="0.8")
    ax.errorbar(radiusne, vsigne, yerr=evsigne, fmt="ro", ecolor="0.8")
    radbin = np.linspace(1, 38, 100)
    # ax.plot(radius, np.polyval(fitall, radius), 'k-')
    # ax.plot(radbin, np.polyval(fitsw, radbin), 'b-')
    # ax.plot(radbin, np.polyval(fitne, radbin), 'r-')
    ax.set_xlabel(r"$R$ (kpc)", size=10)
    ax.set_ylabel(r"$V/\sigma$(R) (km/s)", size=10)
    print figures_dir
    plt.savefig(os.path.join(figures_dir, "vel_over_sigma.png"), dpi=250)


def plot_lambda():
    """ Produces plot for angular momentum parameter lambda as a function of radius. """
    # Selection of PA ranges
    radiusne, lambrne, elambrne = seltabne[:,[0,5,6]].T
    radiussw, lambrsw, elambrsw = seltabsw[:,[0,5,6]].T
    # Fitting
    # fitall = np.polyfit(radius, lambr, 2, w=(1/elambr))
    fitne = np.polyfit(radiusne, lambrne, 3)
    fitsw = np.polyfit(radiussw, lambrsw, 3)
    # Plotting
    plt.style.use("seaborn-paper")
    fs = _large_fig_settings()
    fig = plt.figure(1, figsize=(fs["width"], fs["height"]))
    ax = plt.subplot(1, 1, 1)
    ax.minorticks_on()
    ax.set_xlim(0, 39)
    ax.set_ylim(-0.1, 1.05)
    ax.axvline(x=re, ls="--", c="0.5")
    ax.errorbar(radius, lambr, yerr=elambr, fmt="ko", ecolor="0.8")
    ax.errorbar(radiussw, lambrsw, yerr=elambrsw, fmt="bo", ecolor="0.8")
    ax.errorbar(radiusne, lambrne, yerr=elambrne, fmt="ro", ecolor="0.8")
    radbin = np.linspace(1, 38, 100)
    # ax.plot(radius, np.polyval(fitall, radius), 'k-')
    ax.plot(radbin, np.polyval(fitsw, radbin), 'b-')
    ax.plot(radbin, np.polyval(fitne, radbin), 'r-')
    ax.set_xlabel(r"$R$ (kpc)", size=10)
    ax.set_ylabel(r"$\lambda(R)$", size=10)
    print figures_dir
    plt.savefig(os.path.join(figures_dir, "lambda_r.png"), dpi=250)

def plot_lambdacum():
    """ Produces plot for angular momentum parameter lambda as a function of radius. """
    # Selection of PA ranges
    radiusne, lambcumne = seltabne[:,[0,7]].T
    radiussw, lambcumsw = seltabsw[:,[0,7]].T
    # Fitting
    fitall = np.polyfit(radius, lambcum, 3)
    fitne = np.polyfit(radiusne, lambcumne, 3)
    fitsw = np.polyfit(radiussw, lambcumsw, 3)
    # Plotting
    plt.style.use("seaborn-paper")
    fs = _small_fig_settings()
    fig = plt.figure(1, figsize=(fs["width"], fs["height"]))
    ax = plt.subplot(1, 1, 1)
    ax.minorticks_on()
    ax.set_xlim(0, 39)
    ax.set_ylim(0.0, 0.4)
    ax.axvline(x=re, ls="--", c="0.5")
    ax.plot(radius, lambcum, 'ko')
    # ax.plot(radiusne, lambcumne, 'ro')
    # ax.plot(radiussw, lambcumsw, 'bo')
    radbin = np.linspace(1, 38, 100)
    ax.plot(radius, np.polyval(fitall, radius), 'r-')
    # ax.plot(radbin, np.polyval(fitsw, radbin), 'b-')
    # ax.plot(radbin, np.polyval(fitne, radbin), 'r-')
    ax.set_xlabel(r"$R$ (kpc)", size=12)
    ax.set_ylabel(r"$\lambda(<R)$", size=12)
    plt.subplots_adjust(left=fs["left"], right=fs["right"], bottom=fs["bottom"],
                        top=fs["top"], hspace=fs["hspace"])
    plt.savefig(os.path.join(figures_dir, "lambda_cum_r.png"), dpi=250)

def plot_vsig_lambda():
    """ Produces plot for velocity dispersion as a function of radius. """
    # Selection of PA ranges
    radiussw, vsigsw, evsigsw, lambrsw, elambrsw = seltabsw[:,[0,3,4,5,6]].T
    radiusne, vsigne, evsigne, lambrne, elambrne = seltabne[:,[0,3,4,5,6]].T
    # Fitting v/sigma
    fitallv = np.polyfit(radius, vsig, 3)
    fitnev = np.polyfit(radiusne, vsigne, 3)
    fitswv = np.polyfit(radiussw, vsigsw, 3)
    # Fitting lambda(R)
    fitall = np.polyfit(radius, lambr, 3)
    fitne = np.polyfit(radiusne, lambrne, 3)
    fitsw = np.polyfit(radiussw, lambrsw, 3)
    # Plotting
    plt.style.use("seaborn-paper")
    fs = _large_fig_settings()
    fig = plt.figure(1, figsize=(fs["width"], fs["height"]))
    ax = plt.subplot(2, 1, 1)
    ax.minorticks_on()
    ax.set_xlim(0, 39)
    ax.set_ylim(-0.2, 2.4)
    ax.axvline(x=re, ls="--", c="0.5")
    ax.errorbar(radius, vsig, yerr=evsig, fmt="ko", ecolor="0.8")
    ax.errorbar(radiussw, vsigsw, yerr=evsigsw, fmt="bo", ecolor="0.8")
    ax.errorbar(radiusne, vsigne, yerr=evsigne, fmt="ro", ecolor="0.8")
    radbin = np.linspace(1, 38, 100)
    ax.plot(radbin, np.polyval(fitallv, radbin), 'k-')
    ax.plot(radbin, np.polyval(fitswv, radbin), 'b-')
    ax.plot(radbin, np.polyval(fitnev, radbin), 'r-')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.set_ylabel(r"$V/\sigma$(R) (km/s)", size=10)
    ax = plt.subplot(2, 1, 2)
    ax.minorticks_on()
    ax.set_xlim(0, 39)
    ax.set_ylim(-0.1, 1.05)
    ax.axvline(x=re, ls="--", c="0.5")
    ax.errorbar(radius, lambr, yerr=elambr, fmt="ko", ecolor="0.8")
    ax.errorbar(radiussw, lambrsw, yerr=elambrsw, fmt="bo", ecolor="0.8")
    ax.errorbar(radiusne, lambrne, yerr=elambrne, fmt="ro", ecolor="0.8")
    ax.plot(radbin, np.polyval(fitall, radbin), 'k-')
    ax.plot(radbin, np.polyval(fitsw, radbin), 'b-')
    ax.plot(radbin, np.polyval(fitne, radbin), 'r-')
    ax.set_xlabel(r"$R$ (kpc)", size=10)
    ax.set_ylabel(r"$\lambda(R)$", size=10)
    plt.subplots_adjust(left=fs["left"], right=fs["right"], bottom=fs["bottom"],
                        top=fs["top"], hspace=fs["hspace"])
    plt.savefig(os.path.join(figures_dir, "vsig_lambda_r.png"), dpi=250)

def plot_vsig_lambda2():
    """ Produces plot for velocity dispersion as a function of radius. """
    # Selection of PA ranges
    radiussw, vsigsw, evsigsw, lambrsw, elambrsw, lambcumsw, elambcumsw = seltabsw[:,[0,3,4,5,6,7,8]].T
    radiusne, vsigne, evsigne, lambrne, elambrne, lambcumne, elambcumne = seltabne[:,[0,3,4,5,6,7,8]].T
    # Fitting v/sigma
    fitallv = np.polyfit(radius, vsig, 3)
    # fitnev = np.polyfit(radiusne, vsigne, 3)
    # fitswv = np.polyfit(radiussw, vsigsw, 3)
    # Fitting lambda(R)
    fitall = np.polyfit(radius, lambr, 3)
    # fitne = np.polyfit(radiusne, lambrne, 3)
    # fitsw = np.polyfit(radiussw, lambrsw, 3)
    # Fitting lambda(<R)
    fitallc = np.polyfit(radius, lambcum, 3)
    # fitnec = np.polyfit(radiusne, lambcumne, 3)
    # fitswc = np.polyfit(radiussw, lambcumsw, 3)
    # Plotting
    plt.style.use("seaborn-paper")
    fs = _large_fig_settings()
    fig = plt.figure(1, figsize=(fs["width"], fs["height"]))
    # Top panel
    ax = plt.subplot(3, 1, 1)
    ax.minorticks_on()
    ax.set_xlim(0, 39)
    ax.set_ylim(-0.25, 2.4)
    ax.axvline(x=re, ls="--", c="0.5")
    ax.errorbar(radius, vsig, yerr=evsig, fmt="ko", ecolor="0.8")
    # ax.errorbar(radiussw, vsigsw, yerr=evsigsw, fmt="bo", ecolor="0.8")
    # ax.errorbar(radiusne, vsigne, yerr=evsigne, fmt="ro", ecolor="0.8")
    radbin = np.linspace(1, 38, 100)
    ax.plot(radbin, np.polyval(fitallv, radbin), 'r-')
    # ax.plot(radbin, np.polyval(fitswv, radbin), 'b-')
    # ax.plot(radbin, np.polyval(fitnev, radbin), 'r-')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.set_ylabel(r"$V/\sigma$(R) (km/s)", size=10)
    # Middle panel
    ax = plt.subplot(3, 1, 2)
    ax.minorticks_on()
    ax.set_xlim(0, 39)
    ax.set_ylim(-0.15, 1.08)
    ax.axvline(x=re, ls="--", c="0.5")
    ax.errorbar(radius, lambr, yerr=elambr, fmt="ko", ecolor="0.8")
    # ax.errorbar(radiussw, lambrsw, yerr=elambrsw, fmt="bo", ecolor="0.8")
    # ax.errorbar(radiusne, lambrne, yerr=elambrne, fmt="ro", ecolor="0.8")
    ax.plot(radbin, np.polyval(fitall, radbin), 'r-')
    # ax.plot(radbin, np.polyval(fitsw, radbin), 'b-')
    # ax.plot(radbin, np.polyval(fitne, radbin), 'r-')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.set_ylabel(r"$\lambda(R)$", size=10)
    # Bottom panel
    ax = plt.subplot(3, 1, 3)
    ax.minorticks_on()
    ax.set_xlim(0, 39)
    ax.set_ylim(0.0, 0.42)
    ax.axvline(x=re, ls="--", c="0.5")
    # ax.plot(radius, lambcum, 'ko')
    # ax.plot(radiusne, lambcumne, 'ro')
    # ax.plot(radiussw, lambcumsw, 'bo')
    ax.errorbar(radius, lambcum, yerr=elambcum, fmt="ko", ecolor="0.8")
    # ax.errorbar(radiussw, lambcumsw, yerr=elambcumsw, fmt="bo", ecolor="0.8")
    # ax.errorbar(radiusne, lambcumne, yerr=elambcumne, fmt="ro", ecolor="0.8")
    radbin = np.linspace(1, 38, 100)
    ax.plot(radius, np.polyval(fitallc, radius), 'r-')
    # ax.plot(radbin, np.polyval(fitswc, radbin), 'b-')
    # ax.plot(radbin, np.polyval(fitnec, radbin), 'r-')
    ax.set_xlabel(r"$R$ (kpc)", size=10)
    ax.set_ylabel(r"$\lambda(<R)$", size=10)
    plt.subplots_adjust(left=fs["left"], right=fs["right"], bottom=fs["bottom"],
                        top=fs["top"], hspace=fs["hspace"])
    plt.savefig(os.path.join(figures_dir, "vsig_lambda_cum.png"), dpi=250)

def calc_vsig(vel, evel, sigma, esig, vel0):
    # Calculation of V/sigma
    vsig = (np.abs(vel - vel0)) / sigma
    evsig = vsig * np.sqrt((evel / vel) ** 2 + (esig / sigma) ** 2)
    return vsig, evsig

def calc_lambda(vel, evel, sigma, esig, vel0, radius):
    # Calculation of lambda(R)
    lambr = (radius * np.abs(vel - vel0)) / (radius * np.sqrt((vel - vel0) ** 2 + sigma ** 2))
    elambr = np.sqrt(((sigma ** 2 / ((vel - vel0) ** 2 + sigma ** 2) ** (1.5)) * evel) ** 2
                     + ((-1) * (np.abs(vel - vel0) * sigma) / (((vel - vel0) ** 2 + sigma ** 2) ** (1.5)) * esig) ** 2)
    return lambr, elambr

def calc_lambdacum(vel, evel, sigma, esig, vel0, radius, sn, elambr):
    # Calculation of lambda(<R)
    term1 = sn * radius * np.abs(vel - vel0)
    term1cum = np.cumsum(term1, dtype=float)
    term2 = sn * radius * np.sqrt((vel - vel0)**2 + sigma**2)
    term2cum = np.cumsum(term2, dtype=float)
    lambcum = term1cum/term2cum
    elamb2 = elambr**2
    sumerr = np.cumsum(elamb2, dtype=float)
    edim = np.cumsum((elamb2/elamb2), dtype=float)
    elambcum = np.sqrt(sumerr/edim)
    return lambcum, elambcum


def sel_pane(pane1=90.0, pane2=359.0):
    radiusne = radius[(pa < pane1) | (pa > pane2)]
    sigmane = sigma[(pa < pane1) | (pa > pane2)]
    esigne = esig[(pa < pane1) | (pa > pane2)]
    vsigne = vsig[(pa < pane1) | (pa > pane2)]
    evsigne = evsig[(pa < pane1) | (pa > pane2)]
    lambrne = lambr[(pa < pane1) | (pa > pane2)]
    elambrne = elambr[(pa < pane1) | (pa > pane2)]
    lambcumne = lambcum[(pa < pane1) | (pa > pane2)]
    elambcumne = elambcum[(pa < pane1) | (pa > pane2)]
    seltabne = np.column_stack((radiusne, sigmane, esigne, vsigne, evsigne,
                                lambrne, elambrne, lambcumne, elambcumne))
    return seltabne

def sel_pasw(pasw1=200.0, pasw2=280.0):
    radiussw = radius[(pa > pasw1) & (pa < pasw2)]
    sigmasw = sigma[(pa > pasw1) & (pa < pasw2)]
    esigsw = esig[(pa > pasw1) & (pa < pasw2)]
    vsigsw = vsig[(pa > pasw1) & (pa < pasw2)]
    evsigsw = evsig[(pa > pasw1) & (pa < pasw2)]
    lambrsw = lambr[(pa > pasw1) & (pa < pasw2)]
    elambrsw = elambr[(pa > pasw1) & (pa < pasw2)]
    lambcumsw = lambcum[(pa > pasw1) & (pa < pasw2)]
    elambcumsw = elambcum[(pa > pasw1) & (pa < pasw2)]
    seltabsw = np.column_stack((radiussw, sigmasw, esigsw, vsigsw, evsigsw,
                                lambrsw, elambrsw, lambcumsw, elambcumsw))
    return seltabsw

def _large_fig_settings():
    return {"width":4.0, "height":4.0, "left":0.15, "right":0.99,
            "bottom":0.11, "top":0.99, "hspace":0.1}

def _small_fig_settings():
    return {"width":4.0, "height":2.5, "left":0.13, "right":0.99,
            "bottom":0.12, "top":0.99, "hspace":0.1}

if __name__ == "__main__":
    # Reading tables, calculations and selections
    dataname = os.path.join(results_dir, "n3311_kinematics_sort.dat")
    radius, pa, vel, evel, sigma, esig, sn = np.loadtxt(dataname, usecols=(3, 4, 5, 6, 7, 8, 14), unpack=True)
    vsig, evsig = calc_vsig()
    lambr, elambr = calc_lambda()
    lambcum, elambcum = calc_lambdacum()
    # Selection of PA ranges for NE cone and SW cone
    # pane1 = 90.0 or 153.0; pane2 = 359.0 or 333.0
    # pasw1 = 200.0 or 153.0; pasw2 = 280.0 or 333.0
    seltabne = sel_pane(pane1=90.0, pane2=359.0)
    seltabsw = sel_pasw(pasw1=200.0, pasw2=280.0)
    ###########################################################################
    # Plotting different graphs
    plot_models1()
    # plot_velsigma()
    # plot_lambda()
    # plot_lambdacum()
    # plot_vsig_lambda()
    # plot_vsig_lambda2()
    ###########################################################################