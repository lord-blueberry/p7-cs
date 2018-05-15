
import sys
import os
import numpy as np
import numpy.fft as fourier
import numpy.linalg as nplin
import math
import cmath
from scipy.optimize import minimize

sys.path.append('/home/jon/casa-src/casa/linux/lib/python2.7')
os.environ["CASAPATH"] = "/home/jon/casa-src/casa linux"

from casac import casac as _casac


def MHZToWavelength(x):
    return 299792458.0 / x / 1000000


def inverseFT(vis, u, v, x, y):
    return vis * cmath.exp(-2j*math.pi * (u*x + v*y))


def forwardFT(pixel, u, v, x, y):
    return pixel * cmath.exp(2j*math.pi * (u*x + v*y))


def toComplex(amp, phase):
    x = amp * math.cos(phase)
    y = amp * math.sin(phase)
    return complex(x, y)


image_dimensions = (54, 54)
dirty_api = _casac.image()
dirty_image = dirty_api.newimage(infile="sun_1_vis.image")
dirty_map0 = np.reshape(dirty_image.getchunk(), image_dimensions)[0:54, 0:54]

image_dimensions = (54, 54)
dirty_api2 = _casac.image()
dirty_image2 = dirty_api2.newimage(infile="sun_1_vis_low.image")
dirty_map1 = np.reshape(dirty_image2.getchunk(), image_dimensions)[0:54, 0:54]

dirty_added = dirty_map0 + dirty_map1


pixelIndex = np.zeros(54)
i = 0
for j in range(-27, 27):
    pixelIndex[i] = j
    i += 1
pixelXIndex = pixelIndex * -1 #Because the ZERO of RADEC is bottom right.
pixelYIndex = pixelIndex

MHZ_Base = 1122.0
wave_base = MHZToWavelength(MHZ_Base)

vis0 = toComplex(47.753, math.radians(173.2))
u0 = -112.15 / wave_base
v0 = 576.66 / wave_base

vis1 = toComplex(53.456, math.radians(162.4))
u1 = 71.50 / wave_base
v1 = -368.67 / wave_base

output = np.zeros(image_dimensions)
F_backwards = np.zeros((2, 54*54), dtype=np.complex128)
cell = 8/3600.0
for x in range(0, 54):
    index = x * 54
    for y in range(0, 54):
        x_radians = math.radians((pixelXIndex[x]) * cell)
        y_radians = math.radians((pixelYIndex[y]) * cell)
        F_backwards[0, index + y] = inverseFT(complex(1, 0), u0, v0, x_radians, y_radians)
        F_backwards[1, index + y] = inverseFT(complex(1, 0), u1, v1, x_radians, y_radians)

        tmp0 = inverseFT(vis0, u0, v0, x_radians, y_radians).real
        tmp1 = inverseFT(vis1, u1, v1, x_radians, y_radians).real
        output[x, y] = tmp0 + tmp1
visibilities = np.asarray([vis0, vis1])
output2 = np.real(np.dot(np.asarray([vis0, vis1]), F_backwards))
F = nplin.pinv(F_backwards)
vis_estimates_diff = np.dot(output2.flatten(), F) * 2 - visibilities

axis = 54
image_dimensions= (axis, axis)
pixelIndex = np.zeros(axis)
i = 0
for j in range(-axis/2, axis/2):
    pixelIndex[i] = j
    i += 1
pixelXIndex = pixelIndex * -1 #Because the ZERO of RADEC is bottom right.
pixelYIndex = pixelIndex

inner_try0 = np.zeros(image_dimensions, dtype=np.complex128)
inner_try1 = np.zeros(image_dimensions, dtype=np.complex128)
for x in range(0, image_dimensions[0]):
    for y in range(0, image_dimensions[1]):
        x_radians = math.radians((pixelXIndex[x]) * cell)
        y_radians = math.radians((pixelYIndex[y]) * cell)

        inner_try0[x, y] = forwardFT(1., u0, v0,  x_radians, y_radians)
        inner_try1[x, y] = forwardFT(1., u1, v1, x_radians, y_radians)

vis0_estimate = np.sum(np.multiply(dirty_added, inner_try0)) * 2.0 /inner_try0.size - vis0
vis1_estimate = np.sum(np.multiply(dirty_added, inner_try1)) * 2.0 /inner_try1.size - vis1
print "bla"







#forward transform
sum = complex(0, 0)
for x in range(0, 54):
    for y in range(0, 54):
        x_radians = math.radians((pixelXIndex[x]) * cell)
        y_radians = math.radians((pixelYIndex[y]) * cell)
        sum += forwardFT(dirty_map1[x, y], u1, v1, x_radians, y_radians)


vis1_estimate = sum / output.size
diff = vis1_estimate - vis1

out_api = _casac.image()
out_image = out_api.newimage(infile="sun_1_vis_low.output")
out_reshaped = np.reshape(output, (54, 54, 1, 1))
out_image.putchunk(out_reshaped)
out_image.close()