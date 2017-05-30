# -*- coding: utf-8 -*-
"""

Created on 25/05/2017

@Author: Carlos Eduardo Barbosa

Test of the calculation for the program mh_plot_kinmodels.py from Michael.

"""

import os

import numpy as np
import matplotlib.pyplot as plt

from config import *
from mh_plot_kinmodels import calc_lambda, calc_vsig, calc_lambdacum

def calc_vsig_cb(vel, evel, sigma, esig, vel0):
    # Calculation of V/sigma
    vsig = (np.abs(vel - vel0)) / sigma
    evsig = np.sqrt(evel**2 + (vsig * esig)**2) / sigma
    return vsig, evsig

def calc_lambda_cb(vel, evel, sigma, esig, vel0, radius):
    # Calculation of lambda(R)
    vrel = np.abs(vel - vel0)
    lambr = np.abs(vrel) / np.sqrt(vrel** 2 + sigma**2)
    elambr = np.abs(lambr) * np.sqrt((evel/vrel)**2 + \
                        ((vrel * evel + sigma * esig)/(vrel**2 + sigma**2))**2)
    return lambr, elambr

def calc_lambdacum_cb(vel, evel, sigma, esig, vel0, radius, sn):
    # Calculation of lambda(<R)
    vrel = np.abs(vel - vel0)
    snn = np.where(radius < 25., sn / np.sqrt(3), 1.0 * sn)
    term1 = snn * radius * vrel
    term1cum = np.cumsum(term1, dtype=float)
    term2 = snn * radius * np.sqrt(vrel**2 + sigma**2)
    term2cum = np.cumsum(term2, dtype=float)
    lambcum = term1cum / term2cum
    eterm1 = np.abs(snn * radius * evel)
    eterm1cum = np.sqrt(np.cumsum(eterm1**2, dtype=float))
    eterm2 = np.abs(snn * radius) * np.abs(vrel * evel + sigma * esig) / np.sqrt(vrel**2 + sigma**2)
    eterm2cum = np.sqrt(np.cumsum(eterm2**2, dtype=float))
    elambcum = np.sqrt((eterm1cum / term2cum )**2 + (lambcum * eterm2cum / term2cum)**2)
    return lambcum, elambcum


if __name__ == "__main__":
    table = os.path.join("/home/kadu/Dropbox/hydra1/single2",
                         "n3311_kinematics_sort.dat")
    data = np.loadtxt(table, usecols=(3,5,6,7,8,14))
    data = data[np.argsort(data[:,0])]
    radius, vel, evel, sigma, esig, sn = data.T
    vel0 = np.median(vel)
    # Calculate results from Michael
    vsig_mh, evsig_mh = calc_vsig(vel, evel, sigma, esig, vel0)
    lambr_mh, elambr_mh = calc_lambda(vel, evel, sigma, esig, vel0, radius)
    lambcum_mh, elambcum_mh = calc_lambdacum(vel, evel, sigma, esig, vel0, radius, sn, elambr_mh)
    a_mh = [vsig_mh, evsig_mh, lambr_mh, elambr_mh, lambcum_mh, elambcum_mh]
    # Calculate my results
    vsig, evsig = calc_vsig_cb(vel, evel, sigma, esig, vel0)
    lambr, elambr = calc_lambda_cb(vel, evel, sigma, esig, vel0, radius)
    lambcum, elambcum = calc_lambdacum_cb(vel, evel, sigma, esig, vel0, radius, sn)
    a = [vsig, evsig, lambr, elambr, lambcum, elambcum]
    ax = plt.subplot(231)
    plt.errorbar(radius, vsig, yerr=evsig, fmt='o')
    ax = plt.subplot(232)
    plt.errorbar(radius, lambr, yerr=elambr, fmt='o')
    ax = plt.subplot(233)
    plt.errorbar(radius, lambcum, yerr=elambcum, fmt='o')
    ax = plt.subplot(234)
    plt.errorbar(radius, vsig, yerr=evsig_mh, fmt='o', c="C1")
    ax = plt.subplot(235)
    plt.errorbar(radius, lambr, yerr=elambr_mh, fmt='o', c="C1")
    ax = plt.subplot(236)
    plt.errorbar(radius, lambcum, yerr=elambcum_mh, fmt='o', c="C1")
    plt.show()

