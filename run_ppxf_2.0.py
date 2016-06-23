# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 12:10:54 2014

@author: kadu

Program to run pPXF on hydra I data

Version 1.1: Calculation of S/N is now fixed to a wavelength domain
Version 1.2: Addition of polynomial degrees of pPXF
Version 2.0: Inclusion of MC simulations for the uncertainties
"""
import os
import pickle

import numpy as np
import pyfits as pf
from scipy import ndimage
import matplotlib.pyplot as plt

from ppxf import ppxf
import ppxf_util as util
from config import *
 
def run_ppxf(spectra, velscale):
    """ Run pPXF in a list of spectra"""
    print "Loading templates"
    templates, logLam2, delta = load_templates(velscale)
    for i, spec in enumerate(spectra):
        print "pPXF run of spectrum {0} ({1} of {2})".format(spec, i+1, 
              len(spectra))
        print "Preparing files..."
        # Read a galaxy spectrum and define the wavelength range
        hdu = pf.open(spec)
        spec_lin = hdu[0].data
        h1 = hdu[0].header
        lamRange1 = h1['CRVAL1'] + np.array([0.,h1['CDELT1']*(h1['NAXIS1']-1)])
        # Convolve our spectra to match MILES resolution
        FWHM_dif = np.sqrt(FWHM_tem**2 - FWHM_spec**2)
        sigma = FWHM_dif/2.355/delta # Sigma difference in pixels
        spec_lin = ndimage.gaussian_filter1d(spec_lin,sigma)
        # Rebin to logarithm scale
        galaxy, logLam1, velscale = util.log_rebin(lamRange1, spec_lin, 
                                                   velscale=velscale)
        noise = 0. * galaxy + 0.1
        dv = (logLam2[0]-logLam1[0])*c 
        start, goodPixels = read_setup_file(spec, logLam1)
        print "First interaction of pPXF.."
        pp0 = ppxf(templates, galaxy, noise, velscale, start,
                  goodpixels=goodPixels, plot=False, moments=4,
                  degree=6, mdegree=4, vsyst=dv)
        print "Calculating realistic noise for input..."
        rms0 = galaxy[goodPixels] - pp0.bestfit[goodPixels]
        noise0 = 1.4826 * np.median(np.abs(rms0 - np.median(rms0)))
        noise0 = 0. * galaxy + noise0
        print "Second run of pPXF..."
        pp = ppxf(templates, galaxy, noise0, velscale, pp0.sol,
                  goodpixels=goodPixels, plot=False, moments=4,
                  degree=6, mdegree=4, vsyst=dv, lam=lamRange1)
        print "Finished! Now saving results..."
        with open(spec.replace(".fits", ".pkl"), "w") as f:
            pickle.dump(pp, f)
    return

def wavelength_array(spec):
    """ Produces array for wavelenght of a given array. """
    w0 = pf.getval(spec, "CRVAL1")
    deltaw = pf.getval(spec, "CD1_1")
    pix0 = pf.getval(spec, "CRPIX1")
    npix = pf.getval(spec, "NAXIS1")
    return w0 + deltaw * (np.arange(npix) + 1 - pix0)

def load_templates(velscale):
    """ Load files with stellar library used as templates. """
    current_dir = os.getcwd()
    # Template directory is also set in setup.py
    os.chdir(template_dir)
    miles = [x for x in os.listdir(".") if x.endswith(".fits")]
    miles.sort()
    c = 299792.458
    FWHM_tem = 2.54 # MILES library spectra have a resolution FWHM of 2.54A.
    # Extract the wavelength range and logarithmically rebin one spectrum
    # to the same velocity scale of the SAURON galaxy spectrum, to determine
    # the size needed for the array which will contain the template spectra.
    #
    hdu = pf.open(miles[0])
    ssp = hdu[0].data
    h2 = hdu[0].header
    lamRange2 = h2['CRVAL1'] + np.array([0.,h2['CDELT1']*(h2['NAXIS1']-1)])
    sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp,
                                               velscale=velscale)
    templates = np.empty((sspNew.size,len(miles)))

    for j in range(len(miles)):
        hdu = pf.open(miles[j])
        ssp = hdu[0].data
        sspNew, logLam2, velscale = util.log_rebin(lamRange2, ssp,
                                                   velscale=velscale)
        templates[:,j] = sspNew
    templates /= np.median(templates) # Normalize templates
    os.chdir(current_dir)
    return templates, logLam2, h2['CDELT1']

def read_setup_file(gal, logw):
    w = np.exp(logw)
    filename = gal + ".setup"
    with open(filename) as f:
        f.readline()
        start = f.readline().split()
    start = np.array(start, dtype=float)
    ranges = np.loadtxt(filename, skiprows=5)
    for i, (w1, w2) in enumerate(ranges.reshape((len(ranges)/2, 2))):
        if i == 0:
            good = np.where(np.logical_and(w > w1, w < w2))[0]
        else:
            good = np.hstack((good, np.where(np.logical_and(w > w1,
                                                            w < w2))[0]))
    
    return start, good 

def make_table(spectra, mc=False):
    """ Make table with results. """
    print "Producing summary table..."
    head = ("{0:<30}{1:<14}{2:<14}{3:<14}{4:<14}{5:<14}{6:<14}{7:<14}"
             "{8:<14}{9:<14}{10:<14}{11:<14}{12:<14}{13:<14}\n".format("# FILE",
             "V", "dV", "S", "dS", "h3", "dh3", "h4", "dh4", "chi/DOF",
             "S/N (/ pixel)", "ADEGREE", "MDEGREE", "100*S/N/sigma"))
    results = []
    strmc = "_mc" if mc else ""
    outfile = "ppxf_results{0}.dat".format(strmc)
    for spec in spectra:
        if not os.path.exists(spec.replace(".fits", ".pkl")):
            continue
        pp = pPXF(spec, velscale)
        pp.calc_sn()
        sn = pp.sn
        if mc:
            pp.mc_errors()
        data = [pp.sol[0] + pp.vhelio, pp.error[0],
                pp.sol[1], pp.error[1], pp.sol[2], pp.error[2], pp.sol[3],
                pp.error[3], pp.chi2, sn]
        data = ["{0:12.3f} ".format(x) for x in data]
        line = [spec] + data + ["{0:12}".format(pp.degree),
                                "{0:12}".format(pp.mdegree),
                                "{0:12.3f}".format(100 * sn / pp.sol[1])]
        results.append(line)
    results = np.array(results)
    with open(outfile, "w") as f:
        f.write(head)
        np.savetxt(f, results, fmt="%.30s")
    return

class pPXF():
    """ Class to read pPXF pkl files """
    def __init__(self, spec, velscale):
        """ Load the pkl file from previous ppxf fit and define some attributes.
        """
        with open(spec.replace(".fits", ".pkl")) as f:
            pp = pickle.load(f)
        self.__dict__ = pp.__dict__.copy()
        self.spec = os.path.join(data_dir, spec)
        self.vhelio = pf.getval(self.spec, "VHELIO")
        self.velscale = velscale
        self.w = wavelength_array(os.path.join(data_dir, spec))
        self.flux = pf.getdata(self.spec)
        self.flux_log, self.logw, velscale = util.log_rebin(
                        [self.w[0], self.w[-1]], self.flux, velscale=velscale) 
        self.w_log = np.exp(self.logw)
        return

    def calc_sn(self, w1=5200., w2=5500.):
        idx = np.logical_and(self.w >=w1, self.w <=w2)
        self.res = self.galaxy[idx] - self.bestfit[idx]
        # Using robust method to calculate noise using median deviation
        self.noise = 1.4826 * np.median(np.abs(self.res - np.median(self.res)))
        self.signal = np.sum(self.galaxy[idx]) / len(self.galaxy[idx])
        self.sn = self.signal / self.noise / np.sqrt(0.31)
        return

    def mc_errors(self, nsim=10):
        """ Calculate the errors using MC simulations"""
        errs = np.zeros((nsim, len(self.error)))
        templates, logLam2, delta = load_templates(velscale)
        for i in range(nsim):
            y = self.galaxy + np.random.normal(0, self.noise,
                                                len(self.bestfit))
            noise = np.ones_like(self.galaxy) * self.noise
            sim = ppxf(templates, y, noise, velscale, self.sol,
                       goodpixels=self.goodpixels, plot=False, moments=4,
                       degree=self.degree, mdegree=self.mdegree,
                       vsyst=self.vsyst, lam=self.lam, quiet=True)
            errs[i] = sim.sol
        self.mc_err = np.ma.std(errs, axis=0)
        self.erroor = np.maximum(self.error, self.mc_err)
        return

if __name__ == '__main__':
    # Constants
    c = 299792.458 # Speed of light
    FWHM_tem = 2.54 # MILES library spectra have a resolution FWHM of 2.54A.
    FWHM_spec = 2.1 # FORS2 for Hydra observations has an instrumental
                   # resolution FWHM of 4.2A.
    # Change to data directory according to setup.py program
    # best directory = best .pkl files (a mix of d6m4 and d8m2)
    os.chdir(data_dir)
    # Select spectra to be used in the analysis
    # List of fits files in the folder
    # spectra = [x for x in os.listdir(".") if x.endswith(".pkl")]
    # spectra = [x.replace(".pkl", ".fits") for x in spectra]
    # You can also specify the spectra manually
    # spectra = ["fin_n3311inn2_s38.fits"] # Single spectrum
    spectra = ["fin1_n3311cen1_s14.fits", "fin1_n3311cen2_s45.fits", "fin1_n3311inn2_s39.fits"] # More than one spectrum
    # spectra = ["fin1_n3311cen1_s14.0001.fits", "fin1_n3311cen1_s14.0002.fits", "fin1_n3311cen1_s14.0003.fits", "fin1_n3311cen2_s45.0001.fits", "fin1_n3311cen2_s45.0002.fits", "fin1_n3311cen2_s45.0003.fits", "fin1_n3311inn2_s39.0001.fits", "fin1_n3311inn2_s39.0002.fits", "fin1_n3311inn2_s39.0003.fits", "fin1_n3311inn2_s39.0004.fits", "fin1_n3311inn2_s39.0005.fits"] # More than one
    # Go to the main routine of fitting
    # velscale is defined in the setup.py file, it is used to rebin data
    # run_ppxf(spectra[0:132], velscale) # All spectra
    # run_ppxf(spectra[0:1], velscale) # Single spectrum
    # Make_table produces a table with summary of results and errors
    make_table(spectra, mc=False)
