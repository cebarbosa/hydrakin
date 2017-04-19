# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 10:04:31 2013

@author: kadu

Read headers of spectra to produce tables with info about slitlets.
"""
import os

import pyfits as pf
import numpy as np

from config import *

def field_table_from_header(field, outdir):
    """ Produces info table for a given field. """
    wdir = os.path.join(home, "data/1d")
    # Select one particular spectrum to read the data (this is arbitrary)
    filename = [x for x in os.listdir(wdir) if field in x][0]
    h = pf.getheader(os.path.join(wdir, filename))
    cols = []
    for j in range(101, 200):
        if "ESO INS MOS%s RA" % j in h:
            ra = h["ESO INS MOS%s RA" % j]
            dec = h["ESO INS MOS%s DEC" % j]
            xpos = h["ESO INS MOS%s XPOS" % j]
            ypos = h["ESO INS MOS%s YPOS" % j]
            wid = h["ESO INS MOS%s WID" % j]
            leng = h["ESO INS MOS%s LEN" % j]
            cols.append([ra, dec, xpos, ypos, "n3311" + field, wid, leng])
        else:
            break
    cols = [x for x in cols if x[5] == 1]
    cols = sorted(cols, key=lambda l: l[3], reverse=True)
    outtable = os.path.join(outdir, "n3311{}.txt".format(field))
    with open(outtable, 'w') as g:
        g.write("{0:10}{1:11}{2:12}{3:12}{4:10}{5:12}{6:7}\n".format(
            "# ID", "RA", "DEC", "XPOS", "YPOS", "FIELD", "LENGTH"))
        for k, col in enumerate(cols):
            g.write("{1:4}{0[0]:12.10}{0[1]:12.10}{0[2]:12}{0[3]:12.10}   "
                    "{0[4]:7}{0[6]:6}\n".format(col, k + 1))
    return outtable

def include_types(field, tablename):
    """ Include information about the type of the slit (sky, N3311 or H007)"""
    reftable = os.path.join(tables_dir, "flagged", "n3311{}.txt".format(field))
    refdata = np.loadtxt(reftable, usecols=(5,0,7), dtype="string")
    ids = ["{0[0]}_s{0[1]}".format(x).split("n3311")[1] for x in refdata]
    types = dict(zip(ids, refdata[:,2]))
    with open(tablename) as f:
        lines = [x.strip() for x in f.readlines()]
    header = lines[0] + "{0:>8s}".format("FLAG")
    lines = lines[1:]
    data = np.loadtxt(tablename, usecols=(5, 0), dtype="string")
    dataids = ["{0[0]}_s{0[1]}".format(x).split("n3311")[1] for x in data]
    for i, l in enumerate(lines):
        lines[i] =  l + "{0:>8s}".format(types[dataids[i]])
    with open(tablename, "w") as f:
        f.write(header + "\n")
        f.write("\n".join(lines))
    return

def include_sky_references(tablename):
    sid, xpos, flag = np.loadtxt(tablename, usecols=(0,3,7)).T
    idx0 = np.argwhere(flag == 0).T[0] # Sky slitlets
    idx1 = np.argwhere(flag != 0).T[0] # Non-sky slitlets
    ###########################################################################
    # Calculating references for the non-sky slitlests
    refs = np.zeros_like(idx1)
    for i, idx in enumerate(idx1):
        refs[i] = sid[idx0[np.argsort(np.abs(xpos[idx] - xpos[idx0]))[0]]]
    ###########################################################################
    # Finding which spectra are set for each sky slitlet
    refs0 = []
    for i, idx in enumerate(idx0):
        ref = sid[idx1[np.where(refs == sid[idx])]].astype(int).astype(str
                                                                       ).tolist()
        s = ",".join(ref)
        s = s if s else "None"
        refs0.append(s)
    ###########################################################################
    # Producing new column
    refs0 = np.array(refs0)
    refs = refs.astype(int).astype(str)
    column = np.zeros(len(flag), dtype="S64")
    column[idx0] = refs0
    column[idx1] = refs
    ###########################################################################
    with open(tablename) as f:
        head = f.readline()[:-1]
        lines = [x.strip() for x in f.readlines()]
    for i in range(len(lines)):
        lines[i] = lines[i] +  "     {}".format(column[i])
    with open(tablename, 'w') as g:
        g.write( head + "     {0:>6}\n".format("REF"))
        g.write("\n".join(lines))
    return

def add_extra_slits_n3311(field, outtable):
    if field not in ["cen1", "cen2", "inn2"]:
        return
    with open(outtable) as f:
        head = f.readline()
        data = np.loadtxt(f, dtype=str)
    slits = {"cen1" : ["14", "29"], "cen2" : ["32", "45"], "inn2" : ["39"]}
    names = {"cen1" : [["b", "a", "c"], ["a", "b"]],
             "cen2" : [["a", "b"], ["b", "a", "c"]],
             "inn2" : [["d", "b", "a", "c", "e"]]}
    widths = {"cen1" : [np.array([2.6, 2.8, 3.1]), np.array([3,3])],
              "cen2": [np.array([2.5, 2.5]), np.array([3.3, 3.2, 3.5])],
              "inn2" : [np.array([4.5, 1.5, 1.5, 1.5, 6.5])]}
    theta = np.deg2rad(-40)
    R = np.array([[np.cos(theta), np.sin(theta)],
                  [-np.sin(theta), np.cos(theta)]])
    for slit, ls, ws in zip(slits[field], names[field], widths[field]):
        n = len(ls)
        idx = np.where(data[:, 0] == slit)[0][0]
        line = data[idx]
        ids = np.array([slit + x for x in ls])
        ra = float(line[1])
        dec = float(line[2])
        xpos = np.repeat(line[3], n)
        ypos = np.repeat(line[4], n)
        field = np.repeat(line[5], n)
        length = np.array(ws)
        flag = np.repeat(line[7], n)
        ref = np.repeat(line[8], n)
        #######################################################################
        # Calculating new slitlets
        borders = np.hstack((0, np.cumsum(ws)))
        center = 0.5 * (borders[-1] - borders[0])
        y = borders[:-1] + 0.5 * ws - center
        xy = np.column_stack((np.zeros_like(y), y))
        xy = np.dot(xy, R)
        ras = ra + xy[:, 0] / 3600.
        decs = dec + xy[:, 1] / 3600.
        #######################################################################
        newlines = np.column_stack((ids, ras, decs, xpos, ypos, field,
                                    length, flag, ref))
        data = np.vstack((data, newlines))
        # data[idx,0] = "#" + data[idx,0]
    with open(outtable, "w") as f:
        f.write(head)
        np.savetxt(f, data, fmt="%s")



if __name__ == "__main__":
    fields = ["cen1", "cen2", "inn1", "inn2", "out1", "out2"]
    outdir = os.path.join(tables_dir, "reftables")
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    for field in fields:
        outtable = field_table_from_header(field, outdir)
        include_types(field, outtable)
        include_sky_references(outtable)
        add_extra_slits_n3311(field, outtable)

