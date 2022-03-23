import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval
from lsst.daf import butler as dafButler

from matplotlib import rcParams

rcParams["ytick.labelsize"] = 15
rcParams["xtick.labelsize"] = 15
rcParams["axes.labelsize"] = 20
rcParams["axes.linewidth"] = 2
rcParams["font.size"] = 15
rcParams["axes.titlesize"] = 18

# common functions for AOS analysis


def plotZernikeImage(
    repoDir="/repo/main/",
    collection="u/scichris/Latiss/test",
    instrument="LATISS",
    detector=0,
    titleAdd="",
):
    """Plot the raw Zernikes, postISR image, and
    donut postage stamps.

    Parameters
    ----------
    repoDir : str
        Path to a directory containing butler.yaml file.
        For Latiss tests its /repo/main,  while for
        imgCloseLoop runs it is based on a path to the output
        directory+'phosimData' eg.
        "/project/scichris/aos/rotation_DM-31532/\
         Ns_rotCam_30_c/phosimData/"
    collections: str
        Collection name containing the following dataTypes:
        postISRCCD, zernikeEstimateRaw, donutCatalog

        We assume just one collection containing all data.
        For Latiss tests, collection is u/scichris/Latiss/test

        For ts_phosim runs collections are usually
        f'ts_phosim_90060{iterN}1', where "iterN" is iteration
        number.
    instrument: str
        Name of instrument for which the imgCloseLoop was run,
        eg. LsstCam, LsstComCam, or LATISS
    detector: str or int
        Name of detector. For LATISS it is "0", for
        LsstComCam eg. "R22_S10"; must match the collection name,
        i.e.
    titleAdd: str, optional
        A title to be added to the default main figure title

    Notes
    -------

    Example call to imgCloseLoop:

    python /project/scichris/aos/ts_phosim/bin.src/imgCloseLoop.py \
    --inst comcam --numOfProc 20 --boresightDeg 0.03 -0.02 --rotCam 0 \
    --skyFile /project/scichris/aos/rotation_DM-31532/comCam_grid.txt \
    --output /project/scichris/aos/rotation_DM-31532/Ns_rotCam_0/

    example call to pipetask (running ISR, donut detection, and
    zernike estimation on a pair of auxTel exposures)

    pipetask run -d "exposure IN (2021090800487..2021090800488)  \
    AND instrument='LATISS' AND visit_system=0" \
    -b /repo/main/butler.yaml -i LATISS/raw/all,LATISS/calib  \
    -o u/scichris/Latiss/test \
    -p /project/scichris/aos/ts_wep/tests/testData/\
    pipelineConfigs/testLatissPipeline.yaml --register-dataset-types

    """

    # read in the data from the butler
    butler = dafButler.Butler(repoDir)

    dataId0 = dict(instrument=instrument)
    dataset = next(
        iter(
            butler.registry.queryDatasets(
                datasetType="postISRCCD", dataId=dataId0, collections=[collection]
            )
        )
    )

    expN = dataset.dataId["exposure"]

    # construct a dataId  for postISR
    dataId = {"detector": detector, "instrument": instrument, "exposure": expN}
    print(dataId)
    # read the postISR exposure
    postIsrExp = butler.get("postISRCCD", dataId, collections=[collection])

    # construct a dataId for zernikes and donut catalog:
    # switch exposure to visit
    dataId = {"detector": detector, "instrument": instrument, "visit": expN}
    print(dataId)
    # the raw Zernikes
    zkRaw = butler.get("zernikeEstimateRaw", dataId=dataId, collections=[collection])

    # the donut source catalog
    srcCat = butler.get("donutCatalog", dataId=dataId, collections=[collection])

    # since we queried by detector, sources in that catalog are
    # only for that detector., and its a pandas Df
    exposureName = postIsrExp.getDetector().getName()
    # expCatalog = srcCat.query(f'detector == "{exposureName}"')

    # plot the figure ...
    fig = plt.figure(figsize=(14, 5))

    ##################################
    # left - plot the fit results  ###
    ##################################

    # add_axes([xmin,ymin,dx,dy])
    ax1 = fig.add_axes([0, 0, 0.6, 1])

    for i in range(len(zkRaw)):
        ax1.plot(np.arange(4, 23), 1000 * zkRaw[i], "-d", label=f"donut {i}")

    ax1.set_xlabel(
        "Zernike Number",
    )
    ax1.set_ylabel(
        "Zernike Coefficient [nanometers]",
    )
    ax1.legend(fontsize=14, loc="center left", bbox_to_anchor=[0.65, 0.65])
    ax1.set_xticks(np.arange(4, 23)[::2])
    ax1.grid()

    ax1.set_title(f"{instrument} {collection} {titleAdd}", fontsize=18)

    ##################################
    # right - plot the postISR image #
    ##################################

    ax2 = fig.add_axes([0.6, 0, 0.4, 1])
    exposure_intra = postIsrExp
    zscale = ZScaleInterval()
    data = exposure_intra.image.array
    vmin, vmax = zscale.get_limits(data)

    ax2.imshow(data, origin="lower", vmin=vmin, vmax=vmax)

    nrows = len(srcCat)

    xs = list(srcCat.centroid_x)
    ys = list(srcCat.centroid_y)
    for i in range(nrows):

        x = xs[i]
        y = ys[i]

        # plot the cross marking that the donut was used
        ax2.scatter(x, y, s=200, marker="+", c="m", lw=4)

        # plot the donut number on the plot
        xtext, ytext = x, y
        ytext -= 60
        if xtext + 100 > 4096:
            xtext -= 250
        if len(str(i)) > 1:  # move to the left label thats too long
            # print(i, 'moving')
            xtext -= 340
        else:
            xtext -= 260
        ax2.text(xtext, ytext, f"{i}", fontsize=17, c="white")
    ax2.yaxis.tick_right()
    ax2.set_xlabel("x [px]")
    ax2.set_ylabel("y [px]")
    ax2.yaxis.set_label_position("right")
    ax2.set_title(f"{exposureName}")

    plt.show()

    # plot donuts on a separate figure
    extraFocalStamps = butler.get(
        "donutStampsExtra", dataId=dataId, collections=[collection]
    )
    nDonuts = len(extraFocalStamps)
    ncols = 4
    nrows = nDonuts // ncols
    if nrows * ncols < nDonuts:
        nrows += 1
    fig, axs = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows))
    ax = np.ravel(axs)
    for i in range(nDonuts):
        donut = extraFocalStamps[i]
        ax[i].imshow(donut.stamp_im.image.array, origin="lower")
        ax[i].text(80, 80, f"{i}", fontsize=17, c="white")
    fig.subplots_adjust(hspace=0.35)

    # if there are more axes than donuts,
    # turn off the extra axes
    ncells = nrows * ncols
    if ncells > nDonuts:
        for axis in ax[nDonuts:]:
            axis.axis("off")


def previewExposures(
    yearMonthDay, expStart, expEnd,
    datasetRefOrType="raw",
    collection="LATISS/raw/all",
    instrument="LATISS", detector=0
):
    """Plot auxTel exposures in a given exposure ID range.

    Parameters:
    ----------
    yearMonthDay: str
       The year, month, and day of observations as a string,
       eg. '20210908'.
    expStart: int
       The number of the first exposure in a series to plot, eg. 487.
    expEnd: int
       The number of the final exposure in a series to plot, eg. 490.
    datasetRefOrType: str, optional
        Dataset ref or type for butler, eg. raw or postISRCCD (default: 'raw').
    collection: str, optional
       Collection name containing the auxTel images
       (default: 'LATISS/raw/all').
    instrument: str, optional
        Name of the instrument (default: 'LATISS').
    detector: str or int, optional
        Identifier of a detector (default: 0)

    Notes:
    -------
    It is assumed that all exposures in a chosen range defined by
    range(expStart, expEnd) exist.
    """

    butler = dafButler.Butler("/repo/main/")

    # figure out how many images to plot
    nexp = expEnd - expStart

    # calculate how many cols and rows we need
    if nexp > 3:
        ncol = 3
        nrows = (nexp // ncol) + 1
    else:
        ncol = nexp
        nrows = 1

    zscale = ZScaleInterval()
    # do the plotting
    fig, axs = plt.subplots(nrows, ncol, figsize=(ncol * 4, nrows * 4))
    ax = np.ravel(axs)
    i = 0
    for exp in range(expStart, expEnd):
        exposure = butler.get(
            datasetRefOrType,
            dataId={
                "instrument": instrument,
                "detector": detector,
                "exposure": int(f"{yearMonthDay}00{exp}"),
            },
            collections=[collection],
        )
        data = exposure.image.array
        vmin, vmax = zscale.get_limits(data)
        ax[i].imshow(data, vmin=vmin, vmax=vmax, origin="lower")
        ax[i].set_title(
            f"{yearMonthDay}, exp {exp},\n focusz={np.round(exposure.getMetadata()['FOCUSZ'],3)}"
        )
        i += 1
    fig.subplots_adjust(hspace=0.35)

    # if there are more axes than exposures,
    # turn off the extra axes
    ncells = nrows * ncol
    if ncells > nexp:
        for axis in ax[nexp:]:
            axis.axis("off")
