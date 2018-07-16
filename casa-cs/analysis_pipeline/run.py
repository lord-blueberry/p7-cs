
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
    print("done creating csv images")


def convolve():
    print("nothing")

def createTestImage():
    dimensions = (64, 64)
    image = _casac.image().newimage(infile="./casa_img/sun.flare00_64.64/sun_center_event_64_clean.psf")
    map = np.reshape(image.getchunk(), dimensions)[0:dimensions[0], 0:dimensions[1]]  # weird but this line is needed as is or weird exceptions happen

    test_map = np.zeros(dimensions)
    for i in range(dimensions[0]):
        for j in range(dimensions[1]):
            if i -24 >= 0 and j-13 >=0:
                test_map[i-24, j-13] = map[i, j]

    np.savetxt("test_psf", test_map, delimiter=", ")


def create_centroids():
    inFolder = "/home/jon/Documents/images"
    outFolder = "./centroid_csv"
    algorithms = ["positive_deconv","L1","L2","TV","haar","starlets"]
    dimensions = (64,64)

    for alg in algorithms:
        i = 1
        output_map = np.zeros(dimensions).flatten()
        for timestep in os.listdir(os.path.join(inFolder, alg)):
            for imageName in os.listdir(os.path.join(inFolder, alg, timestep)):
                if imageName.split(".")[-1] == "model":
                    image = _casac.image().newimage(infile=os.path.join(inFolder, alg, timestep, imageName))
                    map = np.reshape(image.getchunk(), dimensions)[0:dimensions[0], 0:dimensions[1]]
                    if map.max() > 0:
                        output_map[np.argmax(map)] = i
                    else:
                        print("img has no max: "+os.path.join(alg,timestep,imageName))
            i = i + 1
        np.savetxt(os.path.join(outFolder,alg+".csv"),np.reshape(output_map,dimensions), delimiter=',')


create_centroids()
#create_csv()

#createTestImage()
#