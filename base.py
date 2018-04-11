
import sys
import os
sys.path.append('/home/jon/casa-src/casa/linux/lib/python2.7')
os.environ["CASAPATH"]="/home/jon/casa-src/casa linux"

from casac import casac as _casac


my_image = _casac.image()
bla = my_image.newimage(infile="/home/jon/raw.img.residual")
myArr = bla.getchunk()
