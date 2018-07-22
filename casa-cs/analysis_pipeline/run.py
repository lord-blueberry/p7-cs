
import sys
import os
import numpy as np

sys.path.append('/home/jon/casa-src/casa/linux/lib/python2.7')
os.environ["CASAPATH"]="/home/jon/casa-src/casa linux"

from casac import casac as _casac

def create_sandbox_csv():
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
    print("done creating csv")

def convert_to_csv(inFolder, outFolder, dimensions):
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)

    for imageName in os.listdir(inFolder):
        imageFile = os.path.join(inFolder, imageName)
        if not (imageFile.endswith(".sumwt") or imageFile.endswith(".csv")):
            print(imageFile)
            image = _casac.image().newimage(infile=imageFile)
            map = np.reshape(image.getchunk(), dimensions)[0:dimensions[0], 0:dimensions[1]]  # weird but this line is needed as is or weird exceptions happen
            np.savetxt(os.path.join(outFolder, imageName+".csv"), map, delimiter=", ")
            image.close()
    print("done creating csv")

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


def create_centroids(inFolder, imgType, outputName):
    outFolder = "./centroid_csv"

    #algorithms = ["positive_deconv","L1","L2","TV","haar","starlets"]
    dimensions_in = (200,200)
    dimensions_out = (100, 100)
    output_map = np.zeros(dimensions_out).flatten()
    for imageFolder in os.listdir(inFolder):
        if imageFolder.split(".")[-1] == imgType:
            folder = os.path.join(inFolder, imageFolder)
            print(folder)
            image = _casac.image().newimage(infile= folder)
            map = np.reshape(image.getchunk(), dimensions_in)[50:150, 50:150]
            if map.max() > 0:
                index = np.argmax(map)
                output_map[index] = output_map[index]+1
            else:
                print("img has no max: "+os.path.join(imageFolder))

    print(np.sum(output_map))
    print(np.max(output_map))
    np.savetxt(os.path.join(outFolder, outputName),np.reshape(output_map,dimensions_out), delimiter=',')


def create_centroids2(inFolder, imgType, outputName):
    outFolder = "./centroid_csv"

    #algorithms = ["positive_deconv","L1","L2","TV","haar","starlets"]
    dimensions = (100,100)
    output_map = np.zeros(dimensions).flatten()
    for imageFolder in os.listdir(inFolder):
        if imageFolder.split(".")[-1] == imgType:
            folder = os.path.join(inFolder, imageFolder)
            print(folder)
            image = _casac.image().newimage(infile= folder)
            map = np.reshape(image.getchunk(), dimensions)[0:dimensions[0], 0:dimensions[1]]
            if map.max() > 0:
                index = np.argmax(map)
                output_map[index] = output_map[index]+1
            else:
                print("img has no max: "+os.path.join(imageFolder))

    print(np.sum(output_map))
    print(np.max(output_map))
    np.savetxt(os.path.join(outFolder, outputName),np.reshape(output_map,dimensions), delimiter=',')

'''
create_centroids("/home/jon/Documents/images/dirty/", "residual","raw.csv")
create_centroids("/home/jon/Documents/images/clean1/", "model","clean1.csv")
create_centroids("/home/jon/Documents/images/clean/", "image","clean.csv")
create_centroids2("/home/jon/Documents/images/peak/", "model","peak.csv")
'''

convert_to_csv("../../results/supernova/casa","../../results/supernova/csv",(128,128))

#create_csv()
#createTestImage()
#