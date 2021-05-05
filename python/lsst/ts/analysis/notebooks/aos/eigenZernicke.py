import numpy as np
from lsst.ts.wep.cwfs.Tool import ZernikeAnnularEval


def createZernikeList(nz, x, y, e):
    '''Create a list of masked zernike images'''
    zimages = []
    for i in range(nz):
        z = np.zeros(nz)
        z[i] = 1.
        zimages.append(ZernikeAnnularEval(z, x, y, e))
    return zimages


def zernikeProjection(w, zernikes, mask):
    '''Calculate zernike coefficients for a given wavefront w'''

    # determine zernike coefficients of wavefront
    nzr = len(zernikes)
    evalues = np.zeros(nzr)
    nMaskPix = len([x for x in mask.flatten() if x == 1])

    for i in range(nzr):
        z1 = zernikes[i]
        dot = np.nansum(w*mask*z1)
        evalues[i] = dot/nMaskPix
    return evalues


def zernikeWeightProjection(w, zernikes, mask):
    '''Calculate zernike coefficients for a wavefront'''

    # determine zernike coefficients of wavefront
    norm = mask.sum()
    nzr = len(zernikes)
    evalues = np.zeros(nzr)
    for i in range(nzr):
        z1 = zernikes[i]
        dot = np.nansum(w*mask*z1)
        evalues[i] = dot/norm
    return evalues


def correlationMatrix(zernikes, mask):
    '''Calculate the correlation matrix giving masked zernikes'''

    # calculate the correlation matrix - assumes masked pixels are nans
    nzr = len(zernikes)
    nMaskPix = len([x for x in mask.flatten() if x == 1])
    correlation = np.zeros((nzr, nzr))

    for i in range(nzr):
        for j in range(nzr):
            z1 = zernikes[i]
            z2 = zernikes[j]
            dot = np.nansum(z1*mask*z2)
            correlation[i][j] = dot
    return correlation/nMaskPix


def correlationWeightMatrix(zernikes, mask):
    '''Calculate the correlation matrix giving masked zernikes'''

    # calculate the correlation matrix - assumes masked pixels are nans
    norm = mask.sum()
    nzr = len(zernikes)
    correlation = np.zeros((nzr, nzr))

    for i in range(nzr):
        for j in range(nzr):
            z1 = zernikes[i]
            z2 = zernikes[j]
            dot = np.nansum(z1*mask*z2)
            correlation[i][j] = dot
    return correlation/norm


def zernikeProjectionMaskCorrection(evalues, correlation):
    '''Correct zernike coefficients for incomplete wavefront'''

    # determine correction and apply to zernike coefficients
    cinv = np.linalg.inv(correlation)
    return cinv@evalues


def reconstructWavefront(evalues, zernikes):
    '''Calculate wavefront from zernike coefficients'''

    nzr = len(evalues)
    reconImage = zernikes[0].copy() * 0.
    for i in range(nzr):
        reconImage += zernikes[i]*evalues[i]
    return reconImage
