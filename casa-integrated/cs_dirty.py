import sys
import os
import numpy as np
import math
import cmath
from scipy.optimize import minimize

sys.path.append('/home/jon/casa-src/casa/linux/lib/python2.7')
os.environ["CASAPATH"] = "/home/jon/casa-src/casa linux"

from casac import casac as _casac

pixel_row = 54
image_dimension = (pixel_row, pixel_row)
folder = "./img/54x54pixels/"
dirty_api = _casac.image()
dirty_image = dirty_api.newimage(infile=folder+"test.residual")
dirty_map = np.reshape(dirty_image.getchunk(), image_dimension)[0:pixel_row, 0:pixel_row]
flat = dirty_map.flatten()

cell = 8/3600.0

#setup forward haar wavelet transform for regularisation


val_lambda = 0.05
def cs(x):
    l2 = np.sum(np.square(np.abs(dirty_map - x)))
    #l1 = np.sum(np.abs(np.dot(haar_transform, x)))
    l1 = 0
    return l2 + val_lambda * l1


x0 = np.zeros(54*54)
result = minimize(cs, x0,method="nelder-mead", options={'xtol': 1e-3, 'disp': True})
np.savetxt("dirty-cs.csv", result.x, delimiter=",")
print "fun: " + str(result.fun)
print "success: " + str(result.success)
