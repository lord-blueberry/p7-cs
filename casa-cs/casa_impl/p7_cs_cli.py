#
# This file was generated using xslt from its XML file
#
# Copyright 2014, Associated Universities Inc., Washington DC
#
import sys
import os
# from casac import *
from casac import casac as _casac
import string
import time
import inspect
import gc
import math
import os
import shutil
from casa_stack_manip import stack_frame_find
from odict import odict
from types import *
from shutil import copyfile

from task_tclean import tclean

import numpy as np
from gurobipy import *

from imagerhelpers.imager_base import PySynthesisImager
from imagerhelpers.input_parameters import ImagerParameters



class p7_cs_cli_:
    __name__ = "p7_cs"
    rkey = None
    i_am_a_casapy_task = None

    @staticmethod
    def copytree(src, dst, symlinks=False, ignore=None):
        if not os.path.exists(dst):
            os.makedirs(dst)
        for item in os.listdir(src):
            s = os.path.join(src, item)
            d = os.path.join(dst, item)
            if os.path.isdir(s):
                copytree(s, d, symlinks, ignore)
            else:
                if not os.path.exists(d) or os.stat(s).st_mtime - os.stat(d).st_mtime > 1:
                    shutil.copy2(s, d)

    @staticmethod
    def read_image(self, imagename, dimensions):
        img = _casac.image()
        dirty = img.newimage(infile=imagename)
        map = np.reshape(dirty.getchunk(), dimensions)[0:dimensions[0], 0:dimensions[1]]
        dirty.close()
        return map

    @staticmethod
    def write_image(self, imagename, content, dimensions):
        img = _casac.image()
        dirty = img.newimage(infile=imagename)
        out_reshaped = np.reshape(content, (dimensions[0], dimensions[1], 1, 1))
        dirty.putchunk(out_reshaped)
        dirty.close()

    @staticmethod
    def model_psf(self, dirty_map, psf_map, psf_threshold, cut_psf=False):
        psf = psf_map.copy()
        psf[np.absolute(psf) < float(psf_threshold)] = 0  # clip negative psf values

        print("psf filled to " + str(float(np.count_nonzero(psf)) / psf.size * 100.0) + "%")

        lo = int(math.ceil(psf.shape[0] / 4))
        hi = psf.shape[0] - int(math.floor(psf.shape[0] / 4))
        if cut_psf:
            psf= psf[lo:hi,lo:hi] #reduce psf size
        psf = np.fliplr(np.flipud(psf))

        #setParam('DistributedMIPJobs',1)
        #setParam('WorkerPool', 'server1062.cs.technik.fhnw.ch:61000')

        model = Model("cs reconstruction")
        model.Params.method = 3 # GRB Primal Simplex == 0   Barrier == 2, concurrent = 3
        psf_sum = model.addVars(dirty_map.shape[0], dirty_map.shape[1], lb=-GRB.INFINITY)
        pixel_flat = []
        pixelArr = []
        for x in range(0, dirty_map.shape[0]):
            row = []
            for y in range(0, dirty_map.shape[1]):
                v = model.addVar()
                row.append(v)
                pixel_flat.append(v)
            pixelArr.append(row)

        XCenter = int(math.floor((psf.shape[0] - 1) / 2))
        YCenter = int(math.floor((psf.shape[1] - 1) / 2))
        print(psf[XCenter, YCenter])
        start_time = time.time()
        for x in range(0, dirty_map.shape[0]):
            psfX0 = -min(x - XCenter, 0)
            psfX1 = min(psf.shape[0] - 1, XCenter + (dirty_map.shape[0] - 1 - x))
            X0 = max(x - XCenter, 0)
            for y in range(0, dirty_map.shape[1]):
                Y0 = max(y - YCenter, 0)
                Y1 = min(y + (psf.shape[1] - YCenter - 1), dirty_map.shape[1] - 1)
                psfY0 = -min(y - YCenter, 0)
                psfY1 = min(psf.shape[1] - 1, YCenter + (dirty_map.shape[1] - 1 - y))

                # print(x, psfX0, psfX1)
                convolution = LinExpr()
                for xp in range(0, psfX1 - psfX0 + 1):
                    psf_cut = psf[xp + psfX0, psfY0:psfY1 + 1]
                    pixel_cut = pixelArr[X0 + xp][Y0:Y1 + 1]
                    convolution.addTerms(psf_cut, pixel_cut)
                model.addConstr(psf_sum[x, y] == dirty_map[x, y] - convolution, "conv")
        elapsed_time = time.time() - start_time
        print("done psf modelling " + str(elapsed_time))

        return model, psf_sum, pixelArr, pixel_flat

    @staticmethod
    def model_peak(self, dirty_map, psf_map, psf_threshold, cut_psf=False):
        model, psf_sum, pixelArr, pixel_flat = self.model_psf(dirty_map, psf_map, psf_threshold, cut_psf)

        model.addSOS(GRB.SOS_TYPE1, pixel_flat)

        objective = QuadExpr()
        for x in range(0, dirty_map.shape[0]):
            for y in range(0, dirty_map.shape[1]):
                # L2[Dirty - X * PSF]
                objective += psf_sum[x, y] * psf_sum[x, y]

        model.setObjective(objective, GRB.MINIMIZE)
        model.optimize()
        objective_val = model.getAttr(GRB.Attr.ObjVal)
        print(objective_val)
        results = np.zeros((dirty_map.shape[0], dirty_map.shape[1]))
        for x in range(0, dirty_map.shape[0]):
            for y in range(0, dirty_map.shape[0]):
                results[x, y] = pixelArr[x][y].x

        return results

    @staticmethod
    def solve_objective_clean(self, dirty_map, psf_map, psf_threshold, cut_psf):
        model, psf_sum, pixelArr, pixel_flat = self.model_psf(dirty_map, psf_map, psf_threshold,cut_psf)

        objective = QuadExpr()
        for x in range(0, dirty_map.shape[0]):
            for y in range(0, dirty_map.shape[1]):
                # L2[Dirty - X * PSF]
                objective += psf_sum[x, y] * psf_sum[x, y]

        model.setObjective(objective, GRB.MINIMIZE)
        model.optimize()
        objective_val = model.getAttr(GRB.Attr.ObjVal)
        print(objective_val)
        results = np.zeros((dirty_map.shape[0], dirty_map.shape[1]))
        for x in range(0, dirty_map.shape[0]):
            for y in range(0, dirty_map.shape[0]):
                results[x, y] = pixelArr[x][y].x

        return results, objective_val

    @staticmethod
    def solve_objective_l1(self, dirty_map, psf_map, psf_threshold, cut_psf, lambda_cs, lambda_estimate=None):
        if lambda_estimate:
            e = np.loadtxt(lambda_estimate[0], delimiter=',')
            model_map = np.loadtxt(lambda_estimate[1], delimiter=',')
            E = np.sum(np.absolute(model_map))
            lambda_cs = e / E
            print("Miller lambda estimation e/E = " + str(lambda_cs))
            print("e: "+str(e))
            print("E: " + str(E))

        model, psf_sum, pixelArr, pixel_flat = self.model_psf(dirty_map, psf_map, psf_threshold)

        abs_pixels = model.addVars(dirty_map.shape[0], dirty_map.shape[1])
        for x in range(0, dirty_map.shape[0]):
            for y in range(dirty_map.shape[1]):
                model.addGenConstrAbs(abs_pixels[x,y], pixelArr[x][y])

        objective = QuadExpr()
        for x in range(0, dirty_map.shape[0]):
            for y in range(0, dirty_map.shape[1]):
                # L2[Dirty - X * PSF]
                objective += psf_sum[x, y] * psf_sum[x, y]

        for x in range(0, dirty_map.shape[0]):
            for y in range(dirty_map.shape[1]):
                objective += lambda_cs * abs_pixels[x,y]

        model.setObjective(objective, GRB.MINIMIZE)
        model.optimize()

        results = np.zeros((dirty_map.shape[0], dirty_map.shape[1]))
        for x in range(0, dirty_map.shape[0]):
            for y in range(0, dirty_map.shape[0]):
                results[x, y] = pixelArr[x][y].x

        return results

    @staticmethod
    def solve_objective_l2(self, dirty_map, psf_map, psf_threshold, cut_psf,lambda_cs,  lambda_estimate=None):
        if lambda_estimate:
            e = np.loadtxt(lambda_estimate[0], delimiter=',')
            model_map = np.loadtxt(lambda_estimate[1], delimiter=',')
            E = np.sum(np.absolute(np.square(model_map)))
            lambda_cs = e / E
            print("Miller lambda estimation " + str(lambda_cs))

        model, psf_sum, pixelArr, pixel_flat = self.model_psf(dirty_map, psf_map, psf_threshold, cut_psf)

        objective = QuadExpr()
        for x in range(0, dirty_map.shape[0]):
            for y in range(0, dirty_map.shape[1]):
                # L2[Dirty - X * PSF]
                objective += psf_sum[x, y] * psf_sum[x, y]

        for x in range(0, dirty_map.shape[0]):
            for y in range(dirty_map.shape[1]):
                objective += lambda_cs * pixelArr[x][y] * pixelArr[x][y]

        model.setObjective(objective, GRB.MINIMIZE)
        model.optimize()


        results = np.zeros((dirty_map.shape[0], dirty_map.shape[1]))
        for x in range(0, dirty_map.shape[0]):
            for y in range(0, dirty_map.shape[0]):
                results[x, y] = pixelArr[x][y].x

        return results

    @staticmethod
    def solve_objective_l1_l2(self, dirty_map, psf_map, psf_threshold, cut_psf, lambda_cs, lambda_estimate=None):
        if lambda_estimate:
            e = np.loadtxt(lambda_estimate[0], delimiter=',')
            model_map = np.loadtxt(lambda_estimate[1], delimiter=',')
            E = np.sum(np.absolute(model_map))
            E = E + np.sum(np.absolute(np.square(model_map)))
            lambda_cs = e / E
            print("Miller lambda estimation e/E = " + str(lambda_cs))
            print("e: " + str(e))
            print("E: " + str(E))

        model, psf_sum, pixelArr, pixel_flat = self.model_psf(dirty_map, psf_map, psf_threshold)

        abs_pixels = model.addVars(dirty_map.shape[0], dirty_map.shape[1])
        for x in range(0, dirty_map.shape[0]):
            for y in range(dirty_map.shape[1]):
                model.addGenConstrAbs(abs_pixels[x, y], pixelArr[x][y])

        objective = QuadExpr()
        for x in range(0, dirty_map.shape[0]):
            for y in range(0, dirty_map.shape[1]):
                # L2[Dirty - X * PSF]
                objective += psf_sum[x, y] * psf_sum[x, y]

        #L2
        for x in range(0, dirty_map.shape[0]):
            for y in range(dirty_map.shape[1]):
                objective += lambda_cs * pixelArr[x][y] * pixelArr[x][y]

        #L1
        for x in range(0, dirty_map.shape[0]):
            for y in range(dirty_map.shape[1]):
                objective += lambda_cs * abs_pixels[x, y]

        model.setObjective(objective, GRB.MINIMIZE)
        model.optimize()

        results = np.zeros((dirty_map.shape[0], dirty_map.shape[1]))
        for x in range(0, dirty_map.shape[0]):
            for y in range(0, dirty_map.shape[0]):
                results[x, y] = pixelArr[x][y].x

        return results

    @staticmethod
    def solve_objective_tv(self, dirty_map, psf_map, psf_threshold, cut_psf,lambda_cs,  lambda_estimate=None):
        if lambda_estimate:
            e = np.loadtxt(lambda_estimate[0], delimiter=',')
            model_map = np. loadtxt(lambda_estimate[1], delimiter=',')
            E = 0
            for x in range(0, dirty_map.shape[0]-1):
                E =  E + np.sum(np.absolute(model_map[x+1] - model_map[x]))
            for y in range(0, dirty_map.shape[1]-1):
                E =  E + np.sum(np.absolute(model_map[:,y+1] - model_map[:,y]))
            lambda_cs = e / E /2
            print("Miller lambda estimation " + str(lambda_cs))

        model, psf_sum, pixelArr, pixel_flat = self.model_psf(dirty_map, psf_map, psf_threshold,cut_psf)

        pixel_x_div = model.addVars(dirty_map.shape[0] - 1, dirty_map.shape[1] - 1, lb=-GRB.INFINITY)
        pixel_x_div_abs = model.addVars(dirty_map.shape[0] - 1, dirty_map.shape[1] - 1, lb=-GRB.INFINITY)
        pixel_y_div = model.addVars(dirty_map.shape[0] - 1, dirty_map.shape[1] - 1, lb=-GRB.INFINITY)
        pixel_y_div_abs = model.addVars(dirty_map.shape[0] - 1, dirty_map.shape[1] - 1, lb=-GRB.INFINITY)
        for x in range(0, dirty_map.shape[0] - 1):
            for y in range(0, dirty_map.shape[1] - 1):
                model.addConstr(pixel_x_div[x, y] == pixelArr[x + 1][y] - pixelArr[x][y])
                model.addGenConstrAbs(pixel_x_div_abs[x, y], pixel_x_div[x, y])

                model.addConstr(pixel_y_div[x, y] == pixelArr[x][y + 1] - pixelArr[x][y])
                model.addGenConstrAbs(pixel_y_div_abs[x, y], pixel_y_div[x, y])

        objective = QuadExpr()
        for x in range(0, dirty_map.shape[0]):
            for y in range(0, dirty_map.shape[1]):
                # L2[Dirty - X * PSF]
                objective += psf_sum[x, y] * psf_sum[x, y]

        for x in range(0, dirty_map.shape[0] - 1):
            for y in range(0, dirty_map.shape[1] - 1):
                objective += lambda_cs * (pixel_x_div_abs[x, y] + pixel_y_div_abs[x, y])

        model.setObjective(objective, GRB.MINIMIZE)
        model.optimize()

        results = np.zeros((dirty_map.shape[0], dirty_map.shape[1]))
        for x in range(0, dirty_map.shape[0]):
            for y in range(0, dirty_map.shape[0]):
                results[x, y] = pixelArr[x][y].x

        return results

    @staticmethod
    def solve_objective_haar(self, dirty_map, psf_map, psf_threshold, cut_psf, lambda_cs,  lambda_estimate=None):
        dimensions = dirty_map.size
        haar2d = np.zeros((dimensions, dimensions))
        norm = np.zeros(dimensions)
        haar2d[0] = 1
        norm[0] = dimensions
        dim = np.zeros(dirty_map.shape)
        def haar_insert(A, norm, x0, y0, x1, y1, k):
            xmid = int((x0 + x1) / 2)
            ymid = int((y0 + y1) / 2)

            dim[x0:xmid, y0:y1] = 1
            dim[xmid:x1, y0:y1] = -1
            haar2d[k] = dim.flatten()
            norm[k] = np.absolute(dim).sum()
            k += 1
            dim[x0:xmid, ymid:y1] = -1
            dim[xmid:x1, ymid:y1] = 1
            haar2d[k] = dim.flatten()
            norm[k] = np.absolute(dim).sum()
            k += 1
            dim[xmid:x1, y0:ymid] = 1
            dim[xmid:x1, ymid:y1] = -1
            haar2d[k] = dim.flatten()
            norm[k] = np.absolute(dim).sum()
            k += 1

            dim[x0:x1, y0:y1] = 0
            if x1 - x0 > 2 and y1 - y0 > 2:
                k = haar_insert(A, norm, x0, y0, xmid, ymid, k)
                k = haar_insert(A, norm, xmid, y0, x1, ymid, k)
                k = haar_insert(A, norm, x0, ymid, xmid, y1, k)
                k = haar_insert(A, norm, xmid, ymid, x1, y1, k)
            return k
        k = 1
        haar_insert(haar2d, norm, 0, 0, dim.shape[0], dim.shape[1], k)
        haar2d = np.transpose(1 / np.sqrt(norm) * np.transpose(haar2d))

        if lambda_estimate:
            e = np.loadtxt(lambda_estimate[0], delimiter=',')
            model_map = np. loadtxt(lambda_estimate[1], delimiter=',')
            E = 0
            flatten = model_map.flatten()
            for x in range(0, dirty_map.size):
                E =  E +np.absolute(np.sum(np.multiply(haar2d[x],flatten)))
            lambda_cs = e / E
            print("Miller lambda estimation " + str(lambda_cs))

        model, psf_sum, pixelArr, pixel_flat = self.model_psf(dirty_map, psf_map, psf_threshold,cut_psf)

        haar_var = model.addVars(dirty_map.size, lb=-GRB.INFINITY)
        haar_var_abs = model.addVars(dirty_map.size, lb=-GRB.INFINITY)
        start_time = time.time()
        for x in range(0, dirty_map.size):
            reg = LinExpr()
            reg.addTerms(haar2d[x], pixel_flat)
            model.addConstr(haar_var[x] == reg)
            model.addGenConstrAbs(haar_var_abs[x], haar_var[x])
        elapsed_time = time.time() - start_time
        print("done regularizer ", elapsed_time)

        objective = QuadExpr()
        for x in range(0, dirty_map.shape[0]):
            offset = x * dirty_map.shape[0]
            for y in range(0, dirty_map.shape[1]):
                # L2[Dirty - X * PSF]
                objective += psf_sum[x, y] * psf_sum[x, y]
                objective += lambda_cs * haar_var_abs[offset + y]

        model.setObjective(objective, GRB.MINIMIZE)
        model.optimize()

        results_haar = np.zeros((dirty_map.shape[0], dirty_map.shape[1]))
        for x in range(0, dirty_map.shape[0]):
            for y in range(0, dirty_map.shape[0]):
                results_haar[x, y] = pixelArr[x][y].x
        return (results_haar)

    @staticmethod
    def solve_objective_starlet(self, dirty_map, psf_map,psf_threshold, cut_psf, lambda_cs, starlet_levels,  lambda_estimate=None):
        casalog = self.__globals__['casalog']
        def calcSpline():
            b3spline = np.asarray([1.0 / 16.0, 1.0 / 4.0, 3.0 / 8.0, 1.0 / 4.0, 1.0 / 16.0])

            row = np.asmatrix([b3spline])
            bla = (np.dot(np.transpose(row), row))
            return (np.dot(np.transpose(row), row))

        def calcConvMatrix(size, kernel, J):
            output = np.zeros((size[0] * size[1], size[0] * size[1]))
            kernel = np.fliplr(np.flipud(kernel))

            disp = 2 ** J
            mid = kernel.shape[0] // 2
            for x in range(0, size[0]):
                offset = x * size[0]
                for y in range(0, size[1]):
                    temp = np.reshape(output[offset + y], size)
                    for i in range(0, kernel.shape[0]):
                        for j in range(0, kernel.shape[1]):
                            xi = ((i - mid) * disp + x)
                            yi = ((j - mid) * disp + y)
                            # print(xi, yi)
                            mx = xi // size[0] % 2
                            xi = xi % size[0]
                            if mx == 1:
                                xi = size[0] - 1 - xi

                            my = yi // size[1] % 2
                            yi = yi % size[1]
                            if my == 1:
                                yi = size[1] - 1 - yi
                            # print(mx, my, xi, yi)
                            temp[xi, yi] += kernel[i, j]
            return (output)

        if lambda_estimate:
            e = np.loadtxt(lambda_estimate[0], delimiter=',')
            model_map = np.loadtxt(lambda_estimate[1], delimiter=',')
            E = 0
            b3spline = calcSpline()
            star_c0 = model_map
            star_cj = 0
            for J in range(0, starlet_levels):
                MJ = calcConvMatrix(dirty_map.shape, b3spline, J)
                star_cj = np.reshape(np.dot(MJ, star_c0.flatten()), model_map.shape)
                star_wj = star_c0 - star_cj
                star_c0 = star_cj
                E = E + np.sum(np.absolute(star_wj))
            E = E + np.sum(np.absolute(star_c0))
            lambda_cs = e / E
            print("Miller lambda estimation "+ str(lambda_cs))

        model, psf_sum, pixelArr, pixel_flat = self.model_psf(dirty_map, psf_map, psf_threshold,cut_psf)

        start_time = time.time()
        star_c = []
        star_w = []
        star_w_abs = []
        b3spline = calcSpline()
        M_previous = 0
        for J in range(0, starlet_levels):
            print("modelling starlet level " + str(J + 1) + " of " + str(starlet_levels))

            star_cJ = []
            star_wJ = []
            star_wJ_abs = []
            for x in range(0, dirty_map.size):
                star_cJ.append(model.addVar(lb=-GRB.INFINITY))
                star_wJ.append(model.addVar(lb=-GRB.INFINITY))
                star_wJ_abs.append(model.addVar())
            star_c.append(star_cJ)
            star_w.append(star_wJ)
            star_w_abs.append(star_wJ_abs)
            star_cminus = 0
            if J == 0:
                star_cminus = pixel_flat
                MJ = calcConvMatrix(dirty_map.shape, b3spline, J)
            else:
                star_cminus = star_c[J - 1]
                MJ = np.dot(M_previous, calcConvMatrix(dirty_map.shape, b3spline, J))
            for x in range(0, dirty_map.size):

                reg = LinExpr()
                reg.addTerms(MJ[x], pixel_flat)
                model.addConstr(star_cJ[x] == reg)
                model.addConstr(star_wJ[x] == star_cminus[x] - star_cJ[x])
                model.addGenConstrAbs(star_wJ_abs[x], star_wJ[x])
            M_previous = MJ
        star_c_abs = []
        for x in range(0, dirty_map.size):
            star_c_abs.append(model.addVar())
            model.addGenConstrAbs(star_c_abs[x], star_c[starlet_levels - 1][x])

        elapsed_time = time.time() - start_time
        print("done starlet modelling ", elapsed_time)

        objective = QuadExpr()
        for x in range(0, dirty_map.shape[0]):
            for y in range(0, dirty_map.shape[1]):
                # L2[Dirty - X * PSF]
                objective += psf_sum[x, y] * psf_sum[x, y]

        for J in range(0, starlet_levels):
            star_wJ_abs = star_w_abs[J]
            for x in range(0, dirty_map.size):
                objective += lambda_cs * star_wJ_abs[x]
        for x in range(0, dirty_map.size):
            objective += lambda_cs * star_c_abs[x]

        model.setObjective(objective, GRB.MINIMIZE)
        model.optimize()

        results_starlet = np.zeros((dirty_map.shape[0], dirty_map.shape[1]))
        for x in range(0, dirty_map.shape[0]):
            for y in range(0, dirty_map.shape[0]):
                results_starlet[x, y] = pixelArr[x][y].x

        return(results_starlet)

    # The existence of the i_am_a_casapy_task attribute allows help()
    # (and other) to treat casapy tasks as a special case.
    def __init__(self):
        self.__bases__ = (p7_cs_cli_,)
        self.__doc__ = self.__call__.__doc__

        self.parameters = {'vis': None, 'selectdata': None, 'field': None, 'spw': None, 'timerange': None,
                           'uvrange': None, 'antenna': None, 'scan': None, 'observation': None, 'intent': None,
                           'datacolumn': None, 'imagename': None, 'imsize': None, 'cell': None, 'phasecenter': None,
                           'stokes': None, 'projection': None, 'startmodel': None, 'specmode': None, 'reffreq': None,
                           'nchan': None, 'start': None, 'width': None, 'outframe': None, 'veltype': None,
                           'restfreq': None, 'interpolation': None, 'gridder': None, 'facets': None, 'chanchunks': None,
                           'wprojplanes': None, 'vptable': None, 'aterm': None, 'psterm': None, 'wbawp': None,
                           'conjbeams': None, 'cfcache': None, 'computepastep': None, 'rotatepastep': None,
                           'pblimit': None, 'normtype': None, 'deconvolver': None, 'scales': None, 'nterms': None,
                           'smallscalebias': None, 'restoration': None, 'restoringbeam': None, 'pbcor': None,
                           'outlierfile': None, 'weighting': None, 'robust': None, 'npixels': None, 'uvtaper': None,
                           'niter': None, 'gain': None, 'threshold': None, 'cycleniter': None, 'cyclefactor': None,
                           'minpsffraction': None, 'maxpsffraction': None, 'interactive': None, 'usemask': None,
                           'mask': None, 'pbmask': None, 'maskthreshold': None, 'maskresolution': None, 'nmask': None,
                           'sidelobethreshold': None, 'noisethreshold': None, 'lownoisethreshold': None,
                           'negativethreshold': None, 'smoothfactor': None, 'minbeamfrac': None, 'cutthreshold': None,
                           'growiterations': None, 'restart': None, 'savemodel': None, 'calcres': None, 'calcpsf': None,
                           'parallel': None,
                           'lambda_cs':None,'cs_alg':None,'psf_threshold':None, 'psf_cutoff':None}

    def result(self, key=None):
        #### and add any that have completed...
        return None

    def __call__(self, vis=None, selectdata=None, field=None, spw=None, timerange=None, uvrange=None, antenna=None,
                 scan=None, observation=None, intent=None, datacolumn=None, imagename=None, imsize=None, cell=None,
                 phasecenter=None, stokes=None, projection=None, startmodel=None, specmode=None, reffreq=None,
                 nchan=None, start=None, width=None, outframe=None, veltype=None, restfreq=None, interpolation=None,
                 gridder=None, facets=None, chanchunks=None, wprojplanes=None, vptable=None, aterm=None, psterm=None,
                 wbawp=None, conjbeams=None, cfcache=None, computepastep=None, rotatepastep=None, pblimit=None,
                 normtype=None, deconvolver=None, scales=None, nterms=None, smallscalebias=None, restoration=None,
                 restoringbeam=None, pbcor=None, outlierfile=None, weighting=None, robust=None, npixels=None,
                 uvtaper=None, niter=None, gain=None, threshold=None, cycleniter=None, cyclefactor=None,
                 minpsffraction=None, maxpsffraction=None, interactive=None, usemask=None, mask=None, pbmask=None,
                 maskthreshold=None, maskresolution=None, nmask=None, sidelobethreshold=None, noisethreshold=None,
                 lownoisethreshold=None, negativethreshold=None, smoothfactor=None, minbeamfrac=None, cutthreshold=None,
                 growiterations=None, restart=None, savemodel=None, calcres=None, calcpsf=None, parallel=None,
                 lambda_cs=None, cs_alg=None, psf_threshold=None, psf_cutoff=None, starlet_levels=None,lambda_estimate=None):

        niter=0 #always zero, we don't want clean to run. we just need the methods to generate the necessary files.

        if not hasattr(self, "__globals__") or self.__globals__ == None:
            self.__globals__ = stack_frame_find()
        # casac = self.__globals__['casac']
        casalog = self.__globals__['casalog']
        casa = self.__globals__['casa']
        # casalog = casac.casac.logsink()
        self.__globals__['__last_task'] = 'p7_cs'
        self.__globals__['taskname'] = 'p7_cs'
        ###
        #self.__globals__['update_params'](func=self.__globals__['taskname'], printtext=False, ipython_globals=self.__globals__)
        ###
        ###
        # Handle globals or user over-ride of arguments
        #
        if type(self.__call__.func_defaults) is NoneType:
            function_signature_defaults = {}
        else:
            function_signature_defaults = dict(
                zip(self.__call__.func_code.co_varnames[1:], self.__call__.func_defaults))
        useLocalDefaults = False

        for item in function_signature_defaults.iteritems():
            key, val = item
            keyVal = eval(key)
            if (keyVal == None):
                # user hasn't set it - use global/default
                pass
            else:
                # user has set it - use over-ride
                if (key != 'self'):
                    useLocalDefaults = True

        myparams = {}
        if useLocalDefaults:
            for item in function_signature_defaults.iteritems():
                key, val = item
                keyVal = eval(key)
                exec ('myparams[key] = keyVal')
                self.parameters[key] = keyVal
                if (keyVal == None):
                    exec ('myparams[key] = ' + key + ' = self.itsdefault(key)')
                    keyVal = eval(key)
                    if (type(keyVal) == dict):
                        if len(keyVal) > 0:
                            exec ('myparams[key] = ' + key + ' = keyVal[len(keyVal)-1][\'value\']')
                        else:
                            exec ('myparams[key] = ' + key + ' = {}')

        else:
            print ''

            myparams['vis'] = vis = self.parameters['vis']
            myparams['selectdata'] = selectdata = self.parameters['selectdata']
            myparams['field'] = field = self.parameters['field']
            myparams['spw'] = spw = self.parameters['spw']
            myparams['timerange'] = timerange = self.parameters['timerange']
            myparams['uvrange'] = uvrange = self.parameters['uvrange']
            myparams['antenna'] = antenna = self.parameters['antenna']
            myparams['scan'] = scan = self.parameters['scan']
            myparams['observation'] = observation = self.parameters['observation']
            myparams['intent'] = intent = self.parameters['intent']
            myparams['datacolumn'] = datacolumn = self.parameters['datacolumn']
            myparams['imagename'] = imagename = self.parameters['imagename']
            myparams['imsize'] = imsize = self.parameters['imsize']
            myparams['cell'] = cell = self.parameters['cell']
            myparams['phasecenter'] = phasecenter = self.parameters['phasecenter']
            myparams['stokes'] = stokes = self.parameters['stokes']
            myparams['projection'] = projection = self.parameters['projection']
            myparams['startmodel'] = startmodel = self.parameters['startmodel']
            myparams['specmode'] = specmode = self.parameters['specmode']
            myparams['reffreq'] = reffreq = self.parameters['reffreq']
            myparams['nchan'] = nchan = self.parameters['nchan']
            myparams['start'] = start = self.parameters['start']
            myparams['width'] = width = self.parameters['width']
            myparams['outframe'] = outframe = self.parameters['outframe']
            myparams['veltype'] = veltype = self.parameters['veltype']
            myparams['restfreq'] = restfreq = self.parameters['restfreq']
            myparams['interpolation'] = interpolation = self.parameters['interpolation']
            myparams['gridder'] = gridder = self.parameters['gridder']
            myparams['facets'] = facets = self.parameters['facets']
            myparams['chanchunks'] = chanchunks = self.parameters['chanchunks']
            myparams['wprojplanes'] = wprojplanes = self.parameters['wprojplanes']
            myparams['vptable'] = vptable = self.parameters['vptable']
            myparams['aterm'] = aterm = self.parameters['aterm']
            myparams['psterm'] = psterm = self.parameters['psterm']
            myparams['wbawp'] = wbawp = self.parameters['wbawp']
            myparams['conjbeams'] = conjbeams = self.parameters['conjbeams']
            myparams['cfcache'] = cfcache = self.parameters['cfcache']
            myparams['computepastep'] = computepastep = self.parameters['computepastep']
            myparams['rotatepastep'] = rotatepastep = self.parameters['rotatepastep']
            myparams['pblimit'] = pblimit = self.parameters['pblimit']
            myparams['normtype'] = normtype = self.parameters['normtype']
            myparams['deconvolver'] = deconvolver = self.parameters['deconvolver']
            myparams['scales'] = scales = self.parameters['scales']
            myparams['nterms'] = nterms = self.parameters['nterms']
            myparams['smallscalebias'] = smallscalebias = self.parameters['smallscalebias']
            myparams['restoration'] = restoration = self.parameters['restoration']
            myparams['restoringbeam'] = restoringbeam = self.parameters['restoringbeam']
            myparams['pbcor'] = pbcor = self.parameters['pbcor']
            myparams['outlierfile'] = outlierfile = self.parameters['outlierfile']
            myparams['weighting'] = weighting = self.parameters['weighting']
            myparams['robust'] = robust = self.parameters['robust']
            myparams['npixels'] = npixels = self.parameters['npixels']
            myparams['uvtaper'] = uvtaper = self.parameters['uvtaper']
            myparams['niter'] = niter = self.parameters['niter']
            myparams['gain'] = gain = self.parameters['gain']
            myparams['threshold'] = threshold = self.parameters['threshold']
            myparams['cycleniter'] = cycleniter = self.parameters['cycleniter']
            myparams['cyclefactor'] = cyclefactor = self.parameters['cyclefactor']
            myparams['minpsffraction'] = minpsffraction = self.parameters['minpsffraction']
            myparams['maxpsffraction'] = maxpsffraction = self.parameters['maxpsffraction']
            myparams['interactive'] = interactive = self.parameters['interactive']
            myparams['usemask'] = usemask = self.parameters['usemask']
            myparams['mask'] = mask = self.parameters['mask']
            myparams['pbmask'] = pbmask = self.parameters['pbmask']
            myparams['maskthreshold'] = maskthreshold = self.parameters['maskthreshold']
            myparams['maskresolution'] = maskresolution = self.parameters['maskresolution']
            myparams['nmask'] = nmask = self.parameters['nmask']
            myparams['sidelobethreshold'] = sidelobethreshold = self.parameters['sidelobethreshold']
            myparams['noisethreshold'] = noisethreshold = self.parameters['noisethreshold']
            myparams['lownoisethreshold'] = lownoisethreshold = self.parameters['lownoisethreshold']
            myparams['negativethreshold'] = negativethreshold = self.parameters['negativethreshold']
            myparams['smoothfactor'] = smoothfactor = self.parameters['smoothfactor']
            myparams['minbeamfrac'] = minbeamfrac = self.parameters['minbeamfrac']
            myparams['cutthreshold'] = cutthreshold = self.parameters['cutthreshold']
            myparams['growiterations'] = growiterations = self.parameters['growiterations']
            myparams['restart'] = restart = self.parameters['restart']
            myparams['savemodel'] = savemodel = self.parameters['savemodel']
            myparams['calcres'] = calcres = self.parameters['calcres']
            myparams['calcpsf'] = calcpsf = self.parameters['calcpsf']
            myparams['parallel'] = parallel = self.parameters['parallel']

        if type(uvtaper) == str: uvtaper = [uvtaper]

        result = None


        paramList = ImagerParameters(
            msname=vis,
            field=field,
            spw=spw,
            timestr=timerange,
            uvdist=uvrange,
            antenna=antenna,
            scan=scan,
            obs=observation,
            state=intent,
            datacolumn=datacolumn,

            ### Image....
            imagename=imagename,
            #### Direction Image Coords
            imsize=imsize,
            cell=cell,
            phasecenter=phasecenter,
            stokes=stokes,
            projection=projection,
            startmodel=startmodel,

            ### Spectral Image Coords
            specmode=specmode,
            reffreq=reffreq,
            nchan=nchan,
            start=start,
            width=width,
            outframe=outframe,
            veltype=veltype,
            restfreq=restfreq,
            sysvel='',  # sysvel,
            sysvelframe='',  # sysvelframe,
            interpolation=interpolation,

            gridder=gridder,
            #        ftmachine=ftmachine,
            facets=facets,
            chanchunks=chanchunks,

            wprojplanes=wprojplanes,

            vptable=vptable,

            ### Gridding....

            aterm=aterm,
            psterm=psterm,
            wbawp=wbawp,
            cfcache=cfcache,
            conjbeams=conjbeams,
            computepastep=computepastep,
            rotatepastep=rotatepastep,

            pblimit=pblimit,
            normtype=normtype,

            outlierfile=outlierfile,
            restart=restart,

            weighting=weighting,
            robust=robust,
            npixels=npixels,
            uvtaper=uvtaper,

            ### Deconvolution
            niter=niter,
            cycleniter=cycleniter,
            loopgain=gain,
            threshold=threshold,
            cyclefactor=cyclefactor,
            minpsffraction=minpsffraction,
            maxpsffraction=maxpsffraction,
            interactive=interactive,

            deconvolver=deconvolver,
            scales=scales,
            nterms=nterms,
            scalebias=smallscalebias,
            restoringbeam=restoringbeam,

            ### new mask params
            usemask=usemask,
            mask=mask,
            pbmask=pbmask,
            maskthreshold=maskthreshold,
            maskresolution=maskresolution,
            nmask=nmask,

            ### automask multithresh params
            sidelobethreshold=sidelobethreshold,
            noisethreshold=noisethreshold,
            lownoisethreshold=lownoisethreshold,
            negativethreshold=negativethreshold,
            smoothfactor=smoothfactor,
            minbeamfrac=minbeamfrac,
            cutthreshold=cutthreshold,
            growiterations=growiterations,

            savemodel=savemodel
        )

        casalog.origin('p7_cs')
        try:
            # if not trec.has_key('tclean') or not casac.casac.utils().verify(mytmp, trec['tclean']) :
            # return False


            scriptstr = ['']
            saveinputs = self.__globals__['saveinputs']
            if type(self.__call__.func_defaults) is NoneType:
                saveinputs = ''
            else:
                saveinputs('tclean', 'tclean.last', myparams, self.__globals__, scriptstr=scriptstr)
            tname = 'p7-cs'
            spaces = ' ' * (18 - len(tname))
            casalog.post('\n##########################################' +
                         '\n##### Begin Task: ' + tname + spaces + ' #####')
            if type(self.__call__.func_defaults) is NoneType:
                casalog.post(scriptstr[0] + '\n', 'INFO')
            else:
                casalog.post(scriptstr[1][1:] + '\n', 'INFO')

            dimensions = (int(imsize[0]), int(imsize[1]))
            print(str(dimensions))

            imager = PySynthesisImager(params=paramList)
            imager.initializeImagers()
            imager.initializeNormalizers()
            imager.setWeighting()

            imager.initializeDeconvolvers()
            imager.initializeIterationControl()

            imager.makePSF()
            imager.makePB()

            imager.runMajorCycle()
            # important magic. otherwise the imager crashes.
            isit = imager.hasConverged()
            imager.runMinorCycle()

            #
            dirty_map = self.read_image(imagename+".residual", dimensions)
            psf_map = self.read_image(imagename + ".psf", dimensions)

            print(imagename)
            if cs_alg =="positive_deconv":
                print("selecting positive deconvolution")
                model_map, objective_val = self.solve_objective_clean(dirty_map, psf_map, psf_threshold, psf_cutoff)
                np.savetxt(imagename+"_objectiveVal.csv", np.asarray(objective_val).reshape(1,), delimiter=',')
                np.savetxt(imagename + "_solution.csv", model_map, delimiter=',')
            elif cs_alg == 'peak':
                print('selecting peak')
                model_map = self.model_peak(dirty_map, psf_map, psf_threshold, psf_cutoff)
            elif cs_alg == "L1":
                print("selecting L1 regularization")
                model_map = self.solve_objective_l1(dirty_map, psf_map, psf_threshold, psf_cutoff, lambda_cs, lambda_estimate)
            elif cs_alg == "L2":
                print("selecting L2 regularization")
                model_map = self.solve_objective_l2(dirty_map, psf_map, psf_threshold, psf_cutoff, lambda_cs, lambda_estimate)
            elif cs_alg == "L1+L2":
                print("selecting L1+L2 mixed regularization")
                model_map = self.solve_objective_l1_l2(dirty_map, psf_map, psf_threshold, psf_cutoff, lambda_cs, lambda_estimate)
            elif cs_alg == "TV":
                print("selecting Total Variation regularization")
                model_map = self.solve_objective_tv(dirty_map, psf_map, psf_threshold, psf_cutoff, lambda_cs, lambda_estimate)
            elif cs_alg == "haar":
                print("selecting haar regularization")
                model_map = self.solve_objective_haar(dirty_map, psf_map, psf_threshold, psf_cutoff, lambda_cs, lambda_estimate)
            else:
                print("selecting starlet regularization")
                if not starlet_levels:
                    print("no starlet levels defined. Setting starlet_levels=1")
                    starlet_levels=1
                model_map = self.solve_objective_starlet(dirty_map, psf_map, psf_threshold, psf_cutoff, lambda_cs, starlet_levels, lambda_estimate)

            self.write_image(imagename+".model", model_map, dimensions)

            imager.runMajorCycle()
            imager.restoreImages()

            imager.deleteTools()

            casalog.post('##### End Task: ' + tname + '  ' + spaces + ' #####' +
                         '\n##########################################')

        except Exception, instance:
            if (self.__globals__.has_key('__rethrow_casa_exceptions') and self.__globals__[
                '__rethrow_casa_exceptions']):
                raise
            else:
                # print '**** Error **** ',instance
                tname = 'p7-cs'
                casalog.post('An error occurred running task ' + tname + '.', 'ERROR')
                raise

        gc.collect()
        return result


    #
    #
    #
    def paramgui(self, useGlobals=True, ipython_globals=None):
        """
        Opens a parameter GUI for this task.  If useGlobals is true, then any relevant global parameter settings are used.
        """
        import paramgui
        if not hasattr(self, "__globals__") or self.__globals__ == None:
            self.__globals__ = stack_frame_find()

        if useGlobals:
            if ipython_globals == None:
                myf = self.__globals__
            else:
                myf = ipython_globals

            paramgui.setGlobals(myf)
        else:
            paramgui.setGlobals({})

        paramgui.runTask('tclean', myf['_ip'])
        paramgui.setGlobals({})

    #
    #
    #
    def defaults(self, param=None, ipython_globals=None, paramvalue=None, subparam=None):
        if not hasattr(self, "__globals__") or self.__globals__ == None:
            self.__globals__ = stack_frame_find()
        if ipython_globals == None:
            myf = self.__globals__
        else:
            myf = ipython_globals

        a = odict()
        a['vis'] = ''
        a['selectdata'] = True
        a['datacolumn'] = 'corrected'
        a['imagename'] = ''
        a['imsize'] = [100]
        a['cell'] = ["1arcsec"]
        a['phasecenter'] = ''
        a['stokes'] = 'I'
        a['projection'] = 'SIN'
        a['startmodel'] = ''
        a['specmode'] = 'mfs'
        a['gridder'] = 'standard'
        a['deconvolver'] = 'hogbom'
        a['restoration'] = True
        a['outlierfile'] = ''
        a['weighting'] = 'natural'
        a['niter'] = 0
        a['usemask'] = 'user'
        a['restart'] = True
        a['savemodel'] = 'none'
        a['calcres'] = True
        a['calcpsf'] = True
        a['parallel'] = False

        a['lambda_cs'] = 0.05
        a['cs_alg'] = "starlet"

        a['selectdata'] = {
            0: odict([{'value': True}, {'field': ""}, {'spw': ""}, {'timerange': ""}, {'uvrange': ""}, {'antenna': ""},
                      {'scan': ""}, {'observation': ""}, {'intent': ""}])}
        a['specmode'] = {
            0: odict([{'value': 'mfs'}, {'reffreq': ""}]),
            1: odict(
                [{'value': 'cube'}, {'nchan': -1}, {'start': ""}, {'width': ""}, {'outframe': ""}, {'veltype': "radio"},
                 {'restfreq': []}, {'interpolation': "linear"}, {'chanchunks': 1}]),
            2: odict([{'value': 'cubedata'}, {'nchan': -1}, {'start': ""}, {'width': ""}, {'veltype': "radio"},
                      {'restfreq': []}, {'interpolation': "linear"}, {'chanchunks': 1}])}
        a['gridder'] = {
            0: odict([{'value': 'standard'}, {'vptable': ""}, {'pblimit': 0.2}]),
            1: odict([{'value': 'widefield'}, {'wprojplanes': 1}, {'facets': 1}, {'vptable': ""}, {'pblimit': 0.2}]),
            2: odict([{'value': 'wproject'}, {'wprojplanes': 1}, {'vptable': ""}, {'pblimit': 0.2}]),
            3: odict([{'value': 'wprojectft'}, {'wprojplanes': 1}, {'vptable': ""}, {'pblimit': 0.2}]),
            4: odict([{'value': 'mosaic'}, {'normtype': "flatnoise"}, {'vptable': ""}, {'pblimit': 0.2},
                      {'conjbeams': False}]),
            5: odict([{'value': 'mosaicft'}, {'normtype': "flatnoise"}, {'vptable': ""}, {'pblimit': 0.2},
                      {'conjbeams': False}]),
            6: odict([{'value': 'ftmosaic'}, {'normtype': "flatnoise"}, {'vptable': ""}, {'pblimit': 0.2}]),
            7: odict([{'value': 'imagemosaic'}, {'wprojplanes': 1}, {'normtype': "flatnoise"}, {'vptable': ""},
                      {'pblimit': 0.2}]),
            8: odict([{'value': 'awproject'}, {'wprojplanes': 1}, {'normtype': "flatnoise"}, {'psterm': False},
                      {'aterm': True}, {'cfcache': ""}, {'computepastep': 360.0}, {'rotatepastep': 360.0},
                      {'wbawp': False}, {'conjbeams': False}, {'pblimit': 0.2}]),
            9: odict([{'value': 'awprojectft'}, {'wprojplanes': 1}, {'normtype': "flatnoise"}, {'psterm': False},
                      {'aterm': True}, {'cfcache': ""}, {'computepastep': 360.0}, {'rotatepastep': 360.0},
                      {'wbawp': False}, {'conjbeams': False}, {'pblimit': 0.2}])}
        a['weighting'] = {
            0: odict([{'value': 'natural'}, {'uvtaper': []}]),
            1: {'value': 'uniform'},
            2: odict([{'value': 'briggs'}, {'robust': 0.5}, {'npixels': 0}, {'uvtaper': []}])}
        a['deconvolver'] = {
            0: {'value': 'hogbom'},
            1: {'value': 'clark'},
            2: odict([{'value': 'multiscale'}, {'scales': []}, {'smallscalebias': 0.6}]),
            3: odict([{'value': 'mtmfs'}, {'scales': []}, {'nterms': 2}]),
            4: {'value': 'aasp'}}
        a['restoration'] = {
            0: odict([{'value': True}, {'restoringbeam': []}, {'pbcor': False}])}
        a['niter'] = {
            0: odict([{'notvalue': 0}, {'gain': 0.1}, {'threshold': 0.0}, {'cycleniter': -1}, {'cyclefactor': 1.0},
                      {'minpsffraction': 0.05}, {'maxpsffraction': 0.8}, {'interactive': False}])}
        a['usemask'] = {
            0: odict([{'value': 'user'}, {'mask': ""}, {'pbmask': 0.0}]),
            1: odict([{'value': 'pb'}, {'pbmask': 0.2}]),
            2: odict([{'value': 'auto-thresh'}, {'pbmask': 0.0}, {'maskthreshold': ""}, {'maskresolution': ""},
                      {'nmask': 0}]),
            3: odict([{'value': 'auto-thresh2'}, {'pbmask': 0.0}, {'maskthreshold': ""}, {'maskresolution': ""},
                      {'nmask': 0}]),
            4: odict(
                [{'value': 'auto-multithresh'}, {'pbmask': 0.2}, {'sidelobethreshold': 3.0}, {'noisethreshold': 5.0},
                 {'lownoisethreshold': 1.5}, {'negativethreshold': 0.0}, {'smoothfactor': 1.0}, {'minbeamfrac': 0.3},
                 {'cutthreshold': 0.01}, {'growiterations': 75}])}

        ### This function sets the default values but also will return the list of
        ### parameters or the default value of a given parameter
        if (param == None):
            myf['__set_default_parameters'](a)
        elif (param == 'paramkeys'):
            return a.keys()
        else:
            if (paramvalue == None and subparam == None):
                if (a.has_key(param)):
                    return a[param]
                else:
                    return self.itsdefault(param)
            else:
                retval = a[param]
                if (type(a[param]) == dict):
                    for k in range(len(a[param])):
                        valornotval = 'value'
                        if (a[param][k].has_key('notvalue')):
                            valornotval = 'notvalue'
                        if ((a[param][k][valornotval]) == paramvalue):
                            retval = a[param][k].copy()
                            retval.pop(valornotval)
                            if (subparam != None):
                                if (retval.has_key(subparam)):
                                    retval = retval[subparam]
                                else:
                                    retval = self.itsdefault(subparam)
                        else:
                            retval = self.itsdefault(subparam)
                return retval

    #
    #
    def check_params(self, param=None, value=None, ipython_globals=None):
        if ipython_globals == None:
            myf = self.__globals__
        else:
            myf = ipython_globals
        #      print 'param:', param, 'value:', value
        try:
            if str(type(value)) != "<type 'instance'>":
                value0 = value
                value = myf['cu'].expandparam(param, value)
                matchtype = False
                if (type(value) == numpy.ndarray):
                    if (type(value) == type(value0)):
                        myf[param] = value.tolist()
                    else:
                        # print 'value:', value, 'value0:', value0
                        # print 'type(value):', type(value), 'type(value0):', type(value0)
                        myf[param] = value0
                        if type(value0) != list:
                            matchtype = True
                else:
                    myf[param] = value
                value = myf['cu'].verifyparam({param: value})
                if matchtype:
                    value = False
        except Exception, instance:
            # ignore the exception and just return it unchecked
            myf[param] = value
        return value

    #
    #
    def description(self, key='tclean', subkey=None):
        desc = {'tclean': 'Radio Interferometric Image Reconstruction',
                'vis': 'Name of input visibility file(s)',
                'selectdata': 'Enable data selection parameters',
                'field': 'field(s) to select',
                'spw': 'spw(s)/channels to select',
                'timerange': 'Range of time to select from data',
                'uvrange': 'Select data within uvrange',
                'antenna': 'Select data based on antenna/baseline',
                'scan': 'Scan number range',
                'observation': 'Observation ID range',
                'intent': 'Scan Intent(s)',
                'datacolumn': 'Data column to image(data,corrected)',
                'imagename': 'Pre-name of output images',
                'imsize': 'Number of pixels',
                'cell': 'Cell size',
                'phasecenter': 'Phase center of the image',
                'stokes': 'Stokes Planes to make',
                'projection': 'Coordinate projection (SIN, HPX)',
                'startmodel': 'Name of starting model image',
                'specmode': 'Spectral definition mode (mfs,cube,cubedata)',
                'reffreq': 'Reference frequency',
                'nchan': 'Number of channels in the output image',
                'start': 'First channel (e.g. start=3,start=\'1.1GHz\',start=\'15343km/s\')',
                'width': 'Channel width (e.g. width=2,width=\'0.1MHz\',width=\'10km/s\')',
                'outframe': 'Spectral reference frame in which to interpret \'start\' and \'width\'',
                'veltype': 'Velocity type (radio, z, ratio, beta, gamma, optical)',
                'restfreq': 'List of rest frequencies',
                'interpolation': 'Spectral interpolation (nearest,linear,cubic)',
                'gridder': 'Gridding options (standard, wproject, widefield, mosaic, awproject)',
                'facets': 'Number of facets on a side',
                'chanchunks': 'Number of channel chunks',
                'wprojplanes': 'Number of distinct w-values for convolution functions',
                'vptable': 'Name of Voltage Pattern table',
                'aterm': 'Use aperture illumination functions during gridding',
                'psterm': 'Use prolate spheroidal during gridding',
                'wbawp': 'Use wideband A-terms',
                'conjbeams': 'Use conjugate frequency for wideband A-terms',
                'cfcache': '>Convolution function cache directory name',
                'computepastep': 'At what parallactic angle interval to recompute AIFs (deg)',
                'rotatepastep': 'At what parallactic angle interval to rotate nearest AIF (deg) ',
                'pblimit': '>PB gain level at which to cut off normalizations ',
                'normtype': 'Normalization type (flatnoise, flatsky)',
                'deconvolver': 'Minor cycle algorithm (hogbom,clark,multiscale,mtmfs,mem,clarkstokes)',
                'scales': 'List of scale sizes (in pixels) for multi-scale algorithms',
                'nterms': 'Number of Taylor coefficients in the spectral model',
                'smallscalebias': 'A bias towards smaller scale sizes',
                'restoration': 'Do restoration steps (or not)',
                'restoringbeam': 'Restoring beam shape to use. Default is the PSF main lobe',
                'pbcor': 'Apply PB correction on the output restored image',
                'outlierfile': 'Name of outlier-field image definitions',
                'weighting': 'Weighting scheme (natural,uniform,briggs)',
                'robust': 'Robustness parameter',
                'npixels': 'Number of pixels to determine uv-cell size (0 : -/+ 3 pixels)',
                'uvtaper': 'uv-taper on outer baselines in uv-plane',
                'niter': 'Maximum number of iterations',
                'gain': 'Loop gain',
                'threshold': 'Stopping threshold ',
                'cycleniter': 'Maximum number of minor-cycle iterations',
                'cyclefactor': 'Scaling on PSF sidelobe level to compute the minor-cycle stopping threshold.',
                'minpsffraction': 'PSF fraction that marks the max depth of cleaning in the minor cycle',
                'maxpsffraction': 'PSF fraction that marks the minimum depth of cleaning in the minor cycle ',
                'interactive': 'Modify masks and parameters at runtime',
                'usemask': 'Type of mask(s) for deconvolution (user, pb, auto-thresh, auto-thresh2, or auto-multithresh)',
                'mask': 'Mask (a list of image name(s) or region file(s) or region string(s) )',
                'pbmask': 'primary beam mask',
                'maskthreshold': 'threshold for automasking (string with unit, e.g. "1.0mJy", sigma,  or fraction of peak ,e.g. 0.1)',
                'maskresolution': 'resolution for automasking (string, e.g. "10arcsec", or a float value as multiplicative factor of the beam)',
                'nmask': 'the maximum number of masks to be added by automasking',
                'sidelobethreshold': 'sidelobethreshold *  the max sidelobe level',
                'noisethreshold': 'noisethreshold * rms in residual image',
                'lownoisethreshold': 'lownoisethreshold * rms in residual image',
                'negativethreshold': 'negativethreshold * rms in residual image',
                'smoothfactor': 'smoothing factor in a unit of the beam',
                'minbeamfrac': 'minimum beam fraction for pruning',
                'cutthreshold': 'threshold to cut the smoothed mask to create a final mask',
                'growiterations': 'number of binary dilation iterations for growing the mask',
                'restart': 'True : Re-use existing images. False : Increment imagename',
                'savemodel': 'Options to save model visibilities (none, virtual, modelcolumn)',
                'calcres': 'Calculate initial residual image',
                'calcpsf': 'Calculate PSF',
                'parallel': 'Run major cycles in parallel',

                'lambda_cs': 'compressed sensing regularization parameter',
                'cs_alg': "setting compressed sensing algorithm",
                'psf_threshold': "psf clip threshold. every pixel with smaller magnitude than psf_threshold will be clipped to 0",
                'psf_cutoff': "cutting off half the psf",
                }

        #
        # Set subfields defaults if needed
        #

        if (desc.has_key(key)):
            return desc[key]

    def itsdefault(self, paramname):
        a = {}
        a['vis'] = ''
        a['selectdata'] = True
        a['field'] = ''
        a['spw'] = ''
        a['timerange'] = ''
        a['uvrange'] = ''
        a['antenna'] = ''
        a['scan'] = ''
        a['observation'] = ''
        a['intent'] = ''
        a['datacolumn'] = 'corrected'
        a['imagename'] = ''
        a['imsize'] = [100]
        a['cell'] = ["1arcsec"]
        a['phasecenter'] = ''
        a['stokes'] = 'I'
        a['projection'] = 'SIN'
        a['startmodel'] = ''
        a['specmode'] = 'mfs'
        a['reffreq'] = ''
        a['nchan'] = -1
        a['start'] = ''
        a['width'] = ''
        a['outframe'] = 'LSRK'
        a['veltype'] = 'radio'
        a['restfreq'] = []
        a['interpolation'] = 'linear'
        a['gridder'] = 'standard'
        a['facets'] = 1
        a['chanchunks'] = 1
        a['wprojplanes'] = 1
        a['vptable'] = ''
        a['aterm'] = True
        a['psterm'] = False
        a['wbawp'] = True
        a['conjbeams'] = False
        a['cfcache'] = ''
        a['computepastep'] = 360.0
        a['rotatepastep'] = 360.0
        a['pblimit'] = 0.2
        a['normtype'] = 'flatnoise'
        a['deconvolver'] = 'hogbom'
        a['scales'] = []
        a['nterms'] = 2
        a['smallscalebias'] = 0.6
        a['restoration'] = True
        a['restoringbeam'] = []
        a['pbcor'] = False
        a['outlierfile'] = ''
        a['weighting'] = 'natural'
        a['robust'] = 0.5
        a['npixels'] = 0
        a['uvtaper'] = ['']
        a['niter'] = 0
        a['gain'] = 0.1
        a['threshold'] = 0.0
        a['cycleniter'] = -1
        a['cyclefactor'] = 1.0
        a['minpsffraction'] = 0.05
        a['maxpsffraction'] = 0.8
        a['interactive'] = False
        a['usemask'] = 'user'
        a['mask'] = ''
        a['pbmask'] = 0.0
        a['maskthreshold'] = ''
        a['maskresolution'] = ''
        a['nmask'] = 0
        a['sidelobethreshold'] = 3.0
        a['noisethreshold'] = 5.0
        a['lownoisethreshold'] = 1.5
        a['negativethreshold'] = 0.0
        a['smoothfactor'] = 1.0
        a['minbeamfrac'] = 0.3
        a['cutthreshold'] = 0.01
        a['growiterations'] = 75
        a['restart'] = True
        a['savemodel'] = 'none'
        a['calcres'] = True
        a['calcpsf'] = True
        a['parallel'] = False

        a['lambda_cs'] = 0.05
        a['cs_alg'] = "starlet"
        a['psf_threshold'] = 0.01
        a['psf_cutoff'] = False

        # a = sys._getframe(len(inspect.stack())-1).f_globals

        if self.parameters['selectdata'] == True:
            a['field'] = ""
            a['spw'] = ""
            a['timerange'] = ""
            a['uvrange'] = ""
            a['antenna'] = ""
            a['scan'] = ""
            a['observation'] = ""
            a['intent'] = ""

        if self.parameters['specmode'] == 'mfs':
            a['reffreq'] = ""

        if self.parameters['specmode'] == 'cube':
            a['nchan'] = -1
            a['start'] = ""
            a['width'] = ""
            a['outframe'] = ""
            a['veltype'] = "radio"
            a['restfreq'] = []
            a['interpolation'] = "linear"
            a['chanchunks'] = 1

        if self.parameters['specmode'] == 'cubedata':
            a['nchan'] = -1
            a['start'] = ""
            a['width'] = ""
            a['veltype'] = "radio"
            a['restfreq'] = []
            a['interpolation'] = "linear"
            a['chanchunks'] = 1

        if self.parameters['gridder'] == 'standard':
            a['vptable'] = ""
            a['pblimit'] = 0.2

        if self.parameters['gridder'] == 'widefield':
            a['wprojplanes'] = 1
            a['facets'] = 1
            a['vptable'] = ""
            a['pblimit'] = 0.2

        if self.parameters['gridder'] == 'wproject':
            a['wprojplanes'] = 1
            a['vptable'] = ""
            a['pblimit'] = 0.2

        if self.parameters['gridder'] == 'wprojectft':
            a['wprojplanes'] = 1
            a['vptable'] = ""
            a['pblimit'] = 0.2

        if self.parameters['gridder'] == 'mosaic':
            a['normtype'] = "flatnoise"
            a['vptable'] = ""
            a['pblimit'] = 0.2
            a['conjbeams'] = False

        if self.parameters['gridder'] == 'mosaicft':
            a['normtype'] = "flatnoise"
            a['vptable'] = ""
            a['pblimit'] = 0.2
            a['conjbeams'] = False

        if self.parameters['gridder'] == 'ftmosaic':
            a['normtype'] = "flatnoise"
            a['vptable'] = ""
            a['pblimit'] = 0.2

        if self.parameters['gridder'] == 'imagemosaic':
            a['wprojplanes'] = 1
            a['normtype'] = "flatnoise"
            a['vptable'] = ""
            a['pblimit'] = 0.2

        if self.parameters['gridder'] == 'awproject':
            a['wprojplanes'] = 1
            a['normtype'] = "flatnoise"
            a['psterm'] = False
            a['aterm'] = True
            a['cfcache'] = ""
            a['computepastep'] = 360.0
            a['rotatepastep'] = 360.0
            a['wbawp'] = False
            a['conjbeams'] = False
            a['pblimit'] = 0.2

        if self.parameters['gridder'] == 'awprojectft':
            a['wprojplanes'] = 1
            a['normtype'] = "flatnoise"
            a['psterm'] = False
            a['aterm'] = True
            a['cfcache'] = ""
            a['computepastep'] = 360.0
            a['rotatepastep'] = 360.0
            a['wbawp'] = False
            a['conjbeams'] = False
            a['pblimit'] = 0.2

        if self.parameters['weighting'] == 'natural':
            a['uvtaper'] = []

        if self.parameters['weighting'] == 'briggs':
            a['robust'] = 0.5
            a['npixels'] = 0
            a['uvtaper'] = []

        if self.parameters['deconvolver'] == 'multiscale':
            a['scales'] = []
            a['smallscalebias'] = 0.6

        if self.parameters['deconvolver'] == 'mtmfs':
            a['scales'] = []
            a['nterms'] = 2

        if self.parameters['restoration'] == True:
            a['restoringbeam'] = []
            a['pbcor'] = False

        if self.parameters['niter'] != 0:
            a['gain'] = 0.1
            a['threshold'] = 0.0
            a['cycleniter'] = -1
            a['cyclefactor'] = 1.0
            a['minpsffraction'] = 0.05
            a['maxpsffraction'] = 0.8
            a['interactive'] = False

        if self.parameters['usemask'] == 'user':
            a['mask'] = ""
            a['pbmask'] = 0.0

        if self.parameters['usemask'] == 'pb':
            a['pbmask'] = 0.2

        if self.parameters['usemask'] == 'auto-thresh':
            a['pbmask'] = 0.0
            a['maskthreshold'] = ""
            a['maskresolution'] = ""
            a['nmask'] = 0

        if self.parameters['usemask'] == 'auto-thresh2':
            a['pbmask'] = 0.0
            a['maskthreshold'] = ""
            a['maskresolution'] = ""
            a['nmask'] = 0

        if self.parameters['usemask'] == 'auto-multithresh':
            a['pbmask'] = 0.2
            a['sidelobethreshold'] = 3.0
            a['noisethreshold'] = 5.0
            a['lownoisethreshold'] = 1.5
            a['negativethreshold'] = 0.0
            a['smoothfactor'] = 1.0
            a['minbeamfrac'] = 0.3
            a['cutthreshold'] = 0.01
            a['growiterations'] = 75

        if a.has_key(paramname):
            return a[paramname]


p7_cs_cli = p7_cs_cli_()