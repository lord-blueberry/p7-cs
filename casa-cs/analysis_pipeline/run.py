
import sys
import os
import numpy as np

sys.path.append('/home/jon/casa-src/casa/linux/lib/python2.7')
os.environ["CASAPATH"]="/home/jon/casa-src/casa linux"

from casac import casac as _casac

def create_csv():
    inFolder = "./casa_img"
    outFolder = "csv_img"
    newNames= {"residual":"dirty.csv", "image":"clean.csv", "psf":"psf.csv", "model":"model.csv"}

    for imageFolder in os.listdir(inFolder):
        outImgFolder = os.path.join(outFolder, imageFolder)
        if not os.path.exists(outImgFolder):
            os.makedirs(outImgFolder)

        strDim = imageFolder.split("_")[-1].split(".")
        dimensions = (int(strDim[0]), int(strDim[1]))
        for imageName in os.listdir(os.path.join(inFolder, imageFolder)):
            newName = newNames[imageName.split(".")[-1]]
            image = _casac.image().newimage(infile=os.path.join(inFolder, imageFolder, imageName))
            map = np.reshape(image.getchunk(), dimensions)[0:dimensions[0], 0:dimensions[1]] #weird but this line is needed as is or weird exceptions happen
            np.savetxt(os.path.join(outImgFolder, newName), map, delimiter=", ")
            image.close()

create_csv()


