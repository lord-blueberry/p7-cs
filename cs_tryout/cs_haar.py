
import sys
import os
import numpy as np
import math
import cmath
from scipy.optimize import minimize

sys.path.append('/home/jon/casa-src/casa/linux/lib/python2.7')
os.environ["CASAPATH"] = "/home/jon/casa-src/casa linux"

from casac import casac as _casac


def MHZToWavelength(x):
    return 299792458.0 / x / 1000000


def inverseFT(vis, u, v, x, y):
    return vis * cmath.exp(-2j*math.pi*(u*x+v*y))


def forwardFT(pixel, u, v, x, y):
    return pixel * cmath.exp(2j*math.pi*(u*x+v*y))


def toComplex(amp, phase):
    x = amp * math.cos(phase)
    y = amp * math.sin(phase)
    return complex(x, y)


pixelIndex = np.zeros(54)
i = 0
for j in range(-27, 27):
    pixelIndex[i] = j
    i += 1
pixelXIndex = pixelIndex * -1
pixelYIndex = pixelIndex

MHZ_Base = 1122.0
wave_base = MHZToWavelength(MHZ_Base)

visibilities = []
u = []
v = []
with open("allvis.sorted.txt") as f:
    for line in f:
        splt = [element for element in line.split(" ") if element != ""]
        amplitude = float(splt[3])
        phase = float(splt[4])
        u_meters = float(splt[9])
        v_meters = float(splt[10])

        visibilities.append(toComplex(amplitude, math.radians(phase)))
        u.append(u_meters / wave_base)
        v.append(v_meters / wave_base)

visibilities = np.array(visibilities)
u = np.array(u)
v = np.array(v)



dirty_api = _casac.image()
dirty_image = dirty_api.newimage(infile="test.residual")
dirty_map = np.reshape(dirty_image.getchunk(), (54, 54))[0:54, 0:54]
flat = dirty_map.flatten()
cell = 8/3600.0
sum = complex(0,0)
for x in range(0, 54):
    for y in range(0, 54):
        x_radians = math.radians((pixelXIndex[x]) * cell)
        y_radians = math.radians((pixelYIndex[y]) * cell)
        tmp = u[275]
        tmp2 = v[275]
        tmp3 = visibilities[275]
        sum += forwardFT(dirty_map[x, y], u[275], v[275], x_radians, y_radians)

bla0 = visibilities[275]
bla1 = sum
diff = visibilities[275] - sum






#setup forward fourier transform matrix
cell = 8/3600.0
F = np.zeros((visibilities.size, 54*54), dtype=np.complex128)
for i in range(0, visibilities.size):
    for x in range(0, 54):
        x_index = x * 54
        for y in range(0, 54):
            x_radians = math.radians((pixelXIndex[x]) * cell)
            y_radians = math.radians((pixelYIndex[y]) * cell)
            F[i, x_index+y] = cmath.exp(2j * math.pi * (u[i] * x_radians + v[i] * y_radians))



#setup forward haar wavelet transform for regularisation
haar_column = np.zeros((54*54))
haar_column[0] = 1/2.0
haar_column[haar_column.size/2] = -1/2.0
haar_transform = np.zeros((54*54, 54*54))
indicator = 0
for i in range(0, haar_column.size):
    haar_transform[i] = haar_column
    indicator += 1
    if indicator == 2:
        indicator = 0
        haar_column = np.roll(haar_column, 1)

haar_transform = np.transpose(haar_transform)



tryout = np.dot(F, flat)

val_lambda = 0.05
def cs(x):
    l2 = np.square(np.abs(visibilities-x*F))
    l1 = np.abs(x * haar_transform)
    return np.sum(l2 + val_lambda * l1)




x0 = np.zeros(54*54)
result = minimize(cs, x0,method="nelder-mead", options={'xtol': 1e-3, 'disp': True})
np.savetxt("my-cs.csv", result.x, delimiter=",")
print "fun: " + str(result.fun)
print "success: " + str(result.success)


