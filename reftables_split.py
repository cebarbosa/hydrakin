# -*- coding: utf-8 -*-
"""

Created on 25/04/16

@author: Carlos Eduardo Barbosa

Split reftables to include sub specs

"""

import os

import numpy as np
import matplotlib.pyplot as plt

from config import *

def split_hcg007(slitdata):
    """ Split data from a given slit in HCC 007 """
    ##########################################################################
    # Center of the sub-slit in relation to the central
    slits = ["cen1_s14", "cen2_s45", "inn2_s39"]
    # Borders defined by Michael's sketch
    borders = [np.array([1., 17 - 5.7, 17 + 5.5, 35.]),
               np.array([1., 21 - 6.2, 21 + 6.2, 40.]),
               np.array([1., 36 - 9., 36 - 3., 36 + 3., 36 + 9., 63.])]
    borders = dict(zip(slits, borders))
    # Identification of the sub slits
    sid= [["b", "a", "c"], ["c", "a", "b"], ["e", "c", "a", "b", "d"]]
    sid = dict(zip(slits, sid))
    ###########################################################################
    # Processing input from table
    id, rac, decc, xpos, ypos, field, leng, flag, ref = slitdata
    mask = field[5:]
    spec = "{0}_s{1}".format(mask, id)
    borders = borders[spec]
    rac = float(rac)
    decc = float(decc)
    xpos = float(xpos)
    ypos = float(ypos)
    leng = float(leng)
    ###########################################################################
    # Calculating new slits
    theta = np.deg2rad(-40)
    borders -= 0.5 * (borders[-1] + borders[0])
    borders = leng * (borders - borders[0]) / (borders[-1] - borders[0])
    sizes = np.diff(borders)
    center = 0.5 * (borders[-1] - borders[0])
    y = borders[:-1] + 0.5 * np.diff(borders) - center
    R = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
    xy = np.column_stack((np.zeros_like(y), y))
    xy = np.dot(xy, R)
    ##########################################################################
    ras = rac  + xy[:,0] / 3600.
    decs = decc + xy[:,1] / 3600.
    lines = []
    for (letter, ra, dec, size) in zip(sid[spec], ras, decs, sizes):
        newid = "{0}{1}".format(spec.split("_")[1][1:], letter)
        newline = "{0} {1} {2} {3} {4:} {5} {6:.1f} {7} " \
                  "{8}".format(newid, ra, dec, xpos, ypos, field,
                               size, flag, ref)
        lines.append(newline.split())
    return np.array(lines)

if __name__ == "__main__":
    reftables_dir = os.path.join(tables_dir, "reftables1")
    output_dir = os.path.join(tables_dir, "reftables2")
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    img_dir = "//home/kadu/Dropbox/hydra1/data/slitimages"
    os.chdir(img_dir)
    masks = ["cen1", "cen2", "inn1", "inn2", "out1", "out2"]
    for mask in masks:
        table = os.path.join(reftables_dir, "n3311{0}.txt".format(mask))
        outtable = os.path.join(output_dir, "n3311{0}.txt".format(mask))
        with open(table) as f:
            head = f.readline()
            data = np.loadtxt(f, dtype=str)
        ######################################################################
        # Append new lines
        idx = np.where(data[:,7]=="3")[0] # Find indices for HCG 007
        if len(idx) > 0:
            newlines = split_hcg007(data[idx[0]])
            # Comment old line
            data[idx[0],0] = "#" + data[idx[0],0]
            data = np.vstack((data, newlines))
        with open(outtable, "w") as f:
            f.write(head)
            np.savetxt(f, data, fmt="%s")