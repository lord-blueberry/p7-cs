
import sys
import os
import numpy as np
import numpy.fft as fourier
import math
import cmath

sys.path.append('/home/jon/casa-src/casa/linux/lib/python2.7')
os.environ["CASAPATH"] = "/home/jon/casa-src/casa linux"

from casac import casac as _casac


def MHZToWavelength(x):
    return 299792458.0 / x / 1000000


def inverseFT(vis, u, v, x, y):
    return vis * cmath.exp(-2j*math.pi*(u*x+v*y))


def toComplex(amp, phase):
    x = amp * math.cos(phase)
    y = amp * math.sin(phase)
    return complex(x, y)


image_dimensions = (54, 54)
dirty_api = _casac.image()
dirty_image = dirty_api.newimage(infile="sun_1_vis_low.image")
dirty_map = np.reshape(dirty_image.getchunk(), image_dimensions)[0:54, 0:54]

pixelIndex = np.zeros(54)
i = 0
for j in range(-26, 28):
    pixelIndex[i] = j
    i += 1
pixelIndex= pixelIndex[::-1]


MHZ_Base = 1122.0
wave_base = MHZToWavelength(MHZ_Base)

vis0 = toComplex(47.753, math.radians(173.2))
u0 = -112.15 / wave_base
v0 = 576.66 / wave_base * -1

vis1 = toComplex(53.456, math.radians(162.4))
u1 = 71.50 / wave_base
v1 = -368.67 / wave_base * -1

output = np.zeros(image_dimensions)
cell = 8/3600.0

for x in range(0, 54):
    for y in range(0, 54):
        x_radians = math.radians((pixelIndex[x]) * cell)
        y_radians = math.radians((pixelIndex[y]) * cell)
        tmp0 = inverseFT(vis0, u0, v0, x_radians, y_radians).real
        tmp1 = inverseFT(vis1, u1, v1, x_radians, y_radians).real
        actual = dirty_map[x, y]
        output[x, y] = tmp0

out_api = _casac.image()
out_image = out_api.newimage(infile="sun_1_vis.output")
out_reshaped = np.reshape(output, (54, 54, 1, 1))
out_image.putchunk(out_reshaped)
out_image.close()