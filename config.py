# -*- coding: utf-8 -*-
"""
Created on 12/09/15

@author: Carlos Eduardo Barbosa

"""

from getpass import getuser

if getuser() == "kadu":
    # Path configurations
    home = "/home/kadu/Dropbox/hydra1"
    template_dir = home + "/templates"
    data_dir = home + "/data/1d/"
    tables_dir = home + "/tables"
    images_dir = home + "/images"
    figures_dir = home + "/figs"
    results_dir = home + "/single2"
else:
    pass

# Set some global constants0
re = 8.4 # Effective radius in kpc, value from Arnaboldi
sn_cut= 10. # Minimum S/N to be used
pa0 = 63. # Photometric position angle
velscale = 30. # Set velocity scale for pPXF related routines

# Constants
c = 299792.458 # Speed of light in km/s
FWHM_tem = 2.54 # MILES library spectra have a resolution FWHM of 2.54A.
FWHM_spec = 2.1 # FORS2 for Hydra observations has an instrumental
               # resolution FWHM of 4.2A.
D = 50.7 # Distance to the center of the Hydra I cluster in Mpc
ra0 = 159.178471651
dec0 = -27.5281283035

