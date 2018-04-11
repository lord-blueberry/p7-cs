
import sys
import os
import numpy as np
import numpy.fft as fourier
import math

sys.path.append('/home/jon/casa-src/casa/linux/lib/python2.7')
os.environ["CASAPATH"]="/home/jon/casa-src/casa linux"

from casac import casac as _casac


dirty_api = _casac.image()
dirty_image = dirty_api.newimage(infile="./img/testImage.residual")
dirty_map = np.reshape(dirty_image.getchunk(), (512, 512))
psf_api = _casac.image()
psf_image = psf_api.newimage(infile="./img/testImage.psf")
psf_map = np.reshape(psf_image.getchunk(), (512, 512))

#spectra helper methods
#!!!!!! Less than and greater than may be reversed!!!
def vis_wv_meyeraux(x):
    y = 35 * x ^ 4 - 84 * x ^ 5 + 70 * x ^ 6 - 20 * x ^ 7
    return y*(x >= 0)*(x <= 1) + (x > 1)

def vis_wv_meyer_wavelet(omega):
    x = abs(omega)
    int1 = ((x > math.Pi/4.) and (x <= math.Pi/2.))
    int2 = ((x > math.Pi/2.) and (x <= math.Pi))
    y = int1 * math.sin(math.Pi/2.*vis_wv_meyeraux(4.*x/math.Pi-1))
    y = y + int2 * math.cos(math.Pi/2*vis_wv_meyeraux(2.*x/math.Pi-1))
    return y

def vis_wv_meyer_scaling(omega):
    x = abs(omega)

    #compute support of Fourier transform of phi.
    int1 = (x < math.Pi / 4.)
    int2 = ((x > math.Pi / 4.) and (x <= math.Pi/ 2.))

    #compute Fourier transform of phi.
    y = int1 + int2 * math.cos(math.Pi  / 2. * vis_wv_meyeraux(4. * x /math.Pi  - 1))

    return y

def vis_wv_linspace(base, limit, n):
    #Provides a row vector V with N linearly spaced elements between BASE and LIMIT;
    v = base + np.arange(n) * (limit - base) / (n - 1)
    return v

def vis_wv_repmat(M0, nc, nr):
    M = M0
    y = 1
    n = max([nc, nr])
    ndims = M.shape.size()
    sm = M.shape

    mx = 1
    my = 1
    if ndims == 1:
        mx = sm[0]
    if ndims > 1:
        mx = sm[0]
        my= sm[1]

    while y < n:
        y = y * 2
        sdmis = M.shape.size()
        sms = M.shape
        smx = 1
        smy = 1
        if ndims == 1:
            smx = sms[0]
        if ndims > 1:
            smx = sms[0]
            smy = sms[1]
        M2 = np.zeros(2 * smx, 2 * smy)
        M2[0:smx - 1, 0:smy - 1] = M
        M2[smx:2 * smx - 1, 0:smy - 1] = M
        M2[0:smx - 1, smy:2 * smy - 1] = M
        M2[smx:2 * smx - 1, smy:2 * smy - 1] = M
        M = M2

    M = M[0:(nc * mx - 1), 0:(nr * my - 1)]
    return M


def vis_wv_fiwt_spectra(imsize, nscales):

    #for better symmetrie each n should be odd
    n_orig = imsize
    n_new = imsize + (1 - (imsize % 2))

    X = 2. ^ (nscales - 1)
    xi_x_init = vis_wv_linspace(0, X, (n_new + 1) / 2)
    #also reverse() of IDL probably interprets 1 as the default axis (probably axis 0 for numpy)
    xi_x_init = [-np.flip(xi_x_init[1:(xi_x_init.size()-1)],1), xi_x_init] #!!!!!! BUG, elements are not the same size. 0 or 1 index begin?
    lx = xi_x_init.size()
    ly = xi_x_init.size()
    xi_x = vis_wv_repmat(xi_x_init,1,ly)
    xi_y = vis_wv_repmat(np.transpose(-np.flip(xi_x_init,1)),lx,1)
    #end bug for sure

    #init
    psi = np.zeros((nscales+1), n_new, n_new)
    psi[0] = vis_wv_meyer_scaling(np.sqrt(xi_x**2 + xi_y**2))
    for j in range(0,nscales-1):
        a = 2.**(-j)
        ax = a * np.sqrt(np.square(xi_x) + np.square(xi_y))
        psi[j + 1] = vis_wv_meyer_wavelet(ax)

    return psi[0:(nscales+1), 0:n_orig, 0:n_orig]




def vis_wv_fiwt(x, psi):
    # foward transform
    nscales = psi.shape.size

    xft = fourier.fftshift(fourier.rfft2(x))
    xft3 = 0 #replicate(complex(0., 0.), nscales, psi.shape[0], psi.shape[1])
    for j in range(0, nscales-1):
        print("bla")
        #xft3[j, *, *] = fft(reform(psi[j, *, *])*xft, / inverse, / center)

    return 0

def vis_wv_fiwt_inverse(x, psi):
    return 0


def vis_wv(vis, imsize, NSCALES=3):
    lamb = 0.05
    psi = vis_wv_fiwt_spectra(imsize[0], NSCALES)
    dc = dirty_map.sum()/dirty_map.size()

    B = dirty_map * (dc / dirty_map.sum())
    P = psf_map/psf_map.sum()

    Btrans= fourier.rfft2(B)

    center=imsize[0]/2
    P = fourier.rfft2(np.roll(P, 1-center))*P.size() #!!!!!
    L = 2 * np.max(np.abs(P) ^ 2)

    #init
    old_total_val = 0
    X_iter = B
    Y = X_iter
    t_new = 1

    for i in range(0, 10):
        X_old = X_iter
        t_old = t_new

        #gradient step
        D = P * fourier.rfft2(Y) - Btrans
        Y = Y -2.0/L * fourier.irfft2(np.conj(P)*D)

        WY = vis_wv_fiwt(Y, psi)

        #soft thresholding
        D = np.abs(WY)- lamb/L
        greaterZero = np.piecewise(D,[D > 0, D <= 0], [1, 0])
        WY = np.sign(np.abs(WY)* (greaterZero * D))

        #the new iterate inverse wavelet transform of WY
        X_iter = vis_wv_fiwt_inverse(WY, psi)

        #flux constraint
        X_iter = X_iter - (np.sum(X_iter)-dc)/ (imsize[0]*imsize[0])

        #updating t and Y
        t_new = (1 + math.sqrt(1. + 4 * t_old ^ 2)) / 2.
        Y = X_iter + ((t_old - 1) / t_new) * (X_iter - X_old)

        #evaluating
        residual = B - fourier.irfft(P*fourier.rfft2(X_iter)) #energy conservation!!!!
        likelyhood = np.sum(np.square(np.abs(residual)))
        sparsity = np.sum(np.abs(vis_wv_fiwt(X_iter, psi)))

        total_val = likelyhood + lamb * sparsity
        old_total_val = total_val

        #stopping criteria
        #if i gt 9 and old_total_val le total_val then break

    return X_iter

imsize = (512, 512)
output = vis_wv(0,imsize,NSCALES=3)

out_api = _casac.image()
out_image = out_api.newimage(infile="./img/testImage.output")
out_reshaped = np.reshape(output, (512, 512, 1, 1))
out_image.putchunk(out_reshaped)
out_image.close()


