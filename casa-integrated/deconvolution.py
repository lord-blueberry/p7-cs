import sys
import os
import numpy as np
import numpy.fft as fourier
import math

sys.path.append('/home/jon/casa-src/casa/linux/lib/python2.7')
os.environ["CASAPATH"]="/home/jon/casa-src/casa linux"

from casac import casac as _casac

#spectra helper methods
def vis_wv_meyeraux(x):
    y = 35 * (x**4) - 84 * (x**5) + 70 * (x**6) - 20 * (x**7)
    return y*(x >= 0)*(x <= 1) + (x > 1)


def vis_wv_meyer_wavelet(omega):
    x = abs(omega)
    int1 = ((x > math.pi/4.) & (x <= math.pi/2.))
    int2 = ((x > math.pi/2.) & (x <= math.pi))
    y = int1 * np.sin(math.pi/2.*vis_wv_meyeraux(4.*x/math.pi-1))
    y = y + int2 * np.cos(math.pi/2*vis_wv_meyeraux(2.*x/math.pi-1))

    return y


def vis_wv_meyer_scaling(omega):
    x = abs(omega)

    #compute support of Fourier transform of phi.
    int1 = (x < math.pi / 4.)
    int2 = (x > math.pi / 4.) & (x <= math.pi/ 2.)

    #compute Fourier transform of phi.
    y = int1 + int2 * np.cos(math.pi / 2. * vis_wv_meyeraux(4. * x /math.pi - 1))
    return y


def vis_wv_linspace(base, limit, n):
    #Provides a row vector V with N linearly spaced elements between BASE and LIMIT;
    v = base + np.arange(n) * (limit - base) / (n - 1)
    return v

#waaaaat? this looks VERY stupid and borderline borked.
def vis_wv_repmat(M0, nc, nr):
    M = M0
    y = 1
    n = max([nc, nr])
    ndims = len(M.shape)
    sm = M.shape

    mx = 1
    my = 1
    if ndims == 1:
        mx = sm[0]
    if ndims > 1:
        mx = sm[0]
        my = sm[1]

    while y < n:
        y = y * 2
        sdmims = len(M.shape)
        sms = M.shape
        smx = 1
        smy = 1
        if sdmims == 1:
            smx = sms[0]
        if sdmims > 1:
            smx = sms[0]
            smy = sms[1]
        M2 = np.zeros((2 * smx, 2 * smy), dtype=np.float64)
        reshaped = np.reshape(M, (smx, smy))
        M2[0:smx, 0:smy] = reshaped
        M2[smx:2 * smx, 0:smy] = reshaped
        M2[0:smx, smy:2 * smy] = reshaped
        M2[smx:2 * smx, smy:2 * smy] = reshaped
        M = M2

    M = M[0:(nc * mx), 0:(nr * my)]
    return M


def vis_wv_fiwt_spectra(imsize, nscales):
    #for better symmetrie each n should be odd
    n_orig = imsize
    n_new = imsize + (1 - (imsize % 2))

    X = 2. ** (nscales - 1)
    xi_x_init = vis_wv_linspace(0, X, (n_new + 1) / 2)
    xi_x_init_flipped = xi_x_init[1:xi_x_init.size]
    xi_x_init_flipped = xi_x_init_flipped[::-1]*-1
    xi_x_init = np.concatenate([xi_x_init_flipped, xi_x_init])
    xi_y_init = xi_x_init[::-1]

    xi_x_init = np.reshape(xi_x_init, (xi_x_init.size, 1))
    xi_y_init = np.reshape(xi_y_init, (xi_x_init.size, 1))
    lx = xi_x_init.size
    ly = xi_y_init.size
    xi_x = vis_wv_repmat(xi_x_init, 1, ly)
    xi_y = vis_wv_repmat(np.transpose(xi_y_init), lx, 1)

    #init
    psi = np.zeros(((nscales+1), n_new, n_new))
    psi[0] = vis_wv_meyer_scaling(np.sqrt(xi_x**2 + xi_y**2))
    for j in range(0, nscales):
        a = 2.**(-j)
        ax = a * np.sqrt(np.square(xi_x) + np.square(xi_y))
        psi[j + 1] = vis_wv_meyer_wavelet(ax)

    return psi[0:(nscales+1), 0:n_orig, 0:n_orig]


def idl_fft(x):
    return fourier.fft2(x)/x.size


def idl_ifft(X):
    return fourier.ifft2(X)*X.size


def vis_wv_fiwt(x, psi):
    # foward transform
    nscales = psi.shape[0]

    xft = fourier.fftshift(idl_fft(x))
    xft3 = np.zeros((nscales, psi.shape[1], psi.shape[2]), dtype=np.complex128)
    for j in range(0, nscales):
        xft3[j] = idl_ifft(fourier.ifftshift(psi[j]*xft))

    return xft3


def vis_wv_fiwt_inverse(x, psi):
    #backwards
    nscales = psi.shape[0]

    r = np.zeros((nscales, psi.shape[1], psi.shape[2]), dtype=np.complex128)
    for j in range(0, nscales):
        r[j] = fourier.fftshift(idl_fft(x[j]))*psi[j] #np.dot instead?

    return idl_ifft(fourier.ifftshift(r.sum(0))) #possible bug: wrong dimension to sum over


def vis_wv(dirty_map, psf_map, imsize, niter, NSCALES=3):
    lamb = 0.05
    psi = vis_wv_fiwt_spectra(imsize[0], NSCALES)
    dc = (dirty_map * (dirty_map > 0)).sum()/dirty_map.size
    dc= 80000

    B = dirty_map * (dc / dirty_map.sum())
    P = psf_map / psf_map.sum()


    Btrans= idl_fft(B)
    center = imsize[0]/2
    trouble = np.transpose(np.roll(np.transpose(P), 1 - center))
    trouble = np.transpose(P)
    P = idl_fft(trouble)*P.size #possible bug
    L = 2 * np.max(abs(P) ** 2)

    #init
    old_total_val = 0
    X_iter = B
    Y = X_iter
    t_new = 1

    for i in range(0, 50):
        X_old = X_iter
        t_old = t_new

        #gradient step
        D = P * idl_fft(Y) - Btrans
        Y = Y# -2.0/L * idl_ifft(np.conj(P)*D)
        WY = vis_wv_fiwt(np.real(Y), psi)

        #soft thresholding
        D = abs(WY) - lamb/L
        WY = np.sign(abs(WY)) * ((D > 0) * D)

        #the new iterate inverse wavelet transform of WY
        X_iter = vis_wv_fiwt_inverse(WY, psi)

        #flux constraint
        X_iter = X_iter - (np.sum(X_iter)-dc) / (imsize[0]*imsize[0])

        #updating t and Yll
        t_new = (1 + math.sqrt(1. + 4 * t_old ** 2)) / 2.
        Y = X_iter + ((t_old - 1) / t_new) * (X_iter - X_old)

        #evaluating
        residual = B - np.real(idl_ifft(P*idl_fft(X_iter))) #energy conservation?
        likelyhood = np.sum(np.square(np.abs(residual)))
        sparsity = np.sum(np.abs(vis_wv_fiwt(X_iter, psi)))

        total_val = likelyhood + lamb * sparsity
        old_total_val = total_val

        #stopping criteria
        #if i gt 9 and old_total_val le total_val then break

    return np.real(X_iter * (X_iter > 0))


folder = "./img/512x512pixels/"
image_dimensions = (512, 512)

dirty_api = _casac.image()
dirty_image = dirty_api.newimage(infile=folder+"test.residual")
dirty_map = np.reshape(dirty_image.getchunk(), image_dimensions)[0:image_dimensions[0], 0:image_dimensions[1]]
psf_api = _casac.image()
psf_image = psf_api.newimage(infile=folder+"test.psf")
psf_map = np.reshape(psf_image.getchunk(), image_dimensions)[0:image_dimensions[0], 0:image_dimensions[1]]

image_dimensions = (53, 53)
output = vis_wv(dirty_map, psf_map, image_dimensions, 50, NSCALES=3)

#output = np.transpose(np.genfromtxt('bla.txt', delimiter=','))
image_dimensions_new=(54, 54)
tmp = np.zeros(image_dimensions_new)
tmp[0:53, 0:53] = output
output = tmp


out_api = _casac.image()
out_image = out_api.newimage(infile=folder+"test.output")
out_reshaped = np.reshape(output, (54, 54, 1, 1))
out_image.putchunk(out_reshaped)
out_image.close()


