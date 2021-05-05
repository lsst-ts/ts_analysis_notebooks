import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rc('image', cmap='jet')


def plotZernikeFits(wavefront, reconImage, trueCoefficients, fitCoefficients, cut=-1.0, title=None):
    '''Plot wavefront images and fits '''

    norm = mpl.colors.Normalize(vmin=-1, vmax=1)

    fig, ax = plt.subplots(1, 3, figsize=(12, 4))
    ax[0].imshow(wavefront)
    ax[1].imshow(wavefront-reconImage)
    ax[1].set_title('{} cut = {}'.format(title, cut))
    ax[2].plot(np.arange(len(trueCoefficients)), trueCoefficients, '-o', label='truth')
    ax[2].plot(np.arange(len(fitCoefficients)), fitCoefficients, 'r*', label='reconstructed')
    ax[2].legend()
    ax[2].grid()
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=None), ax=ax)
