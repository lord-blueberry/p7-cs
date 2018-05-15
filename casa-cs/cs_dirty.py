import sys
import os
import numpy as np
import numpy.fft as fourier
import math
import cmath


sys.path.append('/home/jon/casa-src/casa/linux/lib/python2.7')
os.environ["CASAPATH"] = "/home/jon/casa-src/casa linux"

from casac import casac as _casac

pixel_row = 64
image_dimension = (pixel_row, pixel_row)
folder = "./img/64x64pixels/"
dirty_api = _casac.image()
dirty_image = dirty_api.newimage(infile=folder+"test.residual")
dirty_map = np.reshape(dirty_image.getchunk(), image_dimension)[0:pixel_row, 0:pixel_row]
np.savetxt(folder+"dirty.csv", dirty_map, delimiter=",")

psf_api = _casac.image()
psf_image = psf_api.newimage(infile=folder+"test.psf")
psf_map = np.reshape(psf_image.getchunk(), image_dimension)[0:pixel_row, 0:pixel_row]
np.savetxt(folder+"psf.csv", psf_map, delimiter=",")

clean_api = _casac.image()
clean_image = clean_api.newimage(infile=folder+"test.clean")
clean_map = np.reshape(clean_image.getchunk(), image_dimension)[0:pixel_row, 0:pixel_row]
np.savetxt(folder+"clean.csv", clean_map, delimiter=",")

model_api = _casac.image()
model_image = model_api.newimage(infile=folder+"test.clean")
model_map = np.reshape(model_image.getchunk(), image_dimension)[0:pixel_row, 0:pixel_row]
np.savetxt(folder+"model.csv", model_map, delimiter=",")

