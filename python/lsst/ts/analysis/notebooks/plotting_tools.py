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


def get_zernikes_donuts_from_pickle(zk_results_file):
    """Get raw Zernike fit coefficients and donut stamps
    from pickle file.

    Parameters
    ----------
    zk_results_file: str
        Name of the pickle file with fit results,
        eg. zerDic_2021110400386_extra.npy

    Returns
    -------
    output_zernikes_raw: numpy.ndarray
       Zernike fit coefficients for the donuts.
    donut_stamps: Collection of postage stamps as
            lsst.afw.image.maskedImage.MaskedImage with additional metadata.
    """
    zk_fit = np.load(zk_results_file, allow_pickle=True).item()

    # select only the first element as
    # EstimateZernikesBase stores two
    # identical arrays to match the
    # structure dimension to input dimension
    output_zernikes_raw = zk_fit["outputZernikesRaw"][0]
    donut_stamps = zk_fit["donutStampsExtra"][0]

    return output_zernikes_raw, donut_stamps


def get_zernikes_donuts_postisr_from_butler(repo_dir, collection, instrument, detector):
    """Get raw Zernike fit coefficients, donut stamps,
    and postISR image from butler collection.

    Parameters
    ----------
    repo_dir : str
        Path to a directory containing butler.yaml file.
        For Latiss tests its /repo/main,  while for
        imgCloseLoop runs it is based on a path to the output
        directory+'phosimData' eg.
        "/project/scichris/aos/rotation_DM-31532/\
         Ns_rotCam_30_c/phosimData/"
    collection: str
        Collection name containing the following dataTypes:
        postISRCCD, zernikeEstimateRaw, donutCatalog

        We assume just one collection containing all data.
        For Latiss tests a collection name
        can be "u/scichris/Latiss/test"

        For ts_phosim runs collections are usually
        f'ts_phosim_90060{iterN}1', where "iterN" is iteration
        number.
    instrument: str
        Name of instrument eg. LsstCam, LsstComCam, or LATISS
    detector: str or int
        Name of detector. For LATISS it is 0 , for
        LsstComCam eg. "R22_S10"; must match the collection name
        for imgCloseLoop runs

    Notes
    -----

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

    Returns
    -------
    zernikes_raw: numpy.ndarray
       Zernike fit coefficients for the donuts.
    donut_stamps: Collection of postage stamps as
        lsst.afw.image.maskedImage.MaskedImage.
        with additional metadata.
    exposure: lsst.afw.image.Exposure
        Exposure with the donut image.
    """
    # read in the data from the butler
    butler = dafButler.Butler(repo_dir)

    data_id0 = dict(instrument=instrument)
    dataset = next(
        iter(
            butler.registry.queryDatasets(
                datasetType="postISRCCD", dataId=data_id0, collections=[collection]
            )
        )
    )

    exposure_number = dataset.dataId["exposure"]

    # construct a dataId  for postISR
    data_id = {
        "detector": detector,
        "instrument": instrument,
        "exposure": exposure_number,
    }

    # read the postISR exposure
    exposure = butler.get("postISRCCD", data_id, collections=[collection])

    # construct a dataId for zernikes and donut catalog:
    # switch exposure to visit
    data_id = {"detector": detector, "instrument": instrument, "visit": exposure_number}

    # the raw Zernikes
    zernikes_raw = butler.get(
        "zernikeEstimateRaw", dataId=data_id, collections=[collection]
    )

    donut_stamps = butler.get(
        "donutStampsExtra", dataId=data_id, collections=[collection]
    )

    return zernikes_raw, donut_stamps, exposure


def get_postisr_from_butler(repo_dir, collection, exposure_number):
    """Get the postISR image from butler collection.

    Parameters
    ----------
    repo_dir : str, optional
        Path to a directory containing butler.yaml file.
        For Latiss tests the default is /repo/main
    collection: str
        Collection name containing the postISRCCD.
        Eg. "u/scichris/Latiss/postISRex"
    exposure_number: int
        Exposure number, consisting of year, month, day,
        and sequence number, eg. 2021062300356

    Returns
    -------
    exposure: lsst.afw.image.Exposure
        Exposure with the donut image.
    """
    butler = dafButler.Butler(repo_dir)
    dataset_ref_or_type = "postISRCCD"
    exposure = butler.get(
        dataset_ref_or_type,
        dataId={
            "instrument": "LATISS",
            "detector": 0,
            "exposure": exposure_number,
        },
        collections=[collection],
    )
    return exposure


def plot_raw_zernikes(output_zernikes_raw, ax=None, fig=None, title=""):
    """Plot the raw Zernike coefficients.

    Parameters
    ----------
    zernikes_raw: numpy.ndarray
       Zernike fit coefficients for the donuts.
    ax: axis for plotting, if None it is added to the figure
    fig: figure for plotting
    """
    if fig is None:
        fig = plt.figure()

    if ax is None:
        ax = fig.add_axes([0, 0, 0.6, 1])

    for i in range(len(output_zernikes_raw)):

        ax.plot(
            np.arange(4, 23), 1000 * output_zernikes_raw[i], "-d", label=f"donut {i}"
        )

    ax.set_xlabel(
        "Zernike Number",
    )
    ax.set_ylabel(
        "Zernike Coefficient [nanometers]",
    )
    ax.legend(fontsize=14, loc="center left", bbox_to_anchor=[0.55, 0.65], ncol=2)
    ax.set_xticks(np.arange(4, 23)[::2])
    ax.grid()

    ax.set_title(title, fontsize=18)

    return


def plot_donut_locations(
    exposure,
    donut_stamps,
    ax=None,
    fig=None,
):
    """Plot the donut locations on the postISR image.

    Parameters
    ----------
    exposure: lsst.afw.image.Exposure
        Exposure (postISR) with the donut images.
    donut_stamps: Collection of postage stamps as
        lsst.afw.image.maskedImage.MaskedImage.
        with additional metadata.
    ax: axis for plotting, if None it is added to the figure.
    fig: figure for plotting , if None a new figure is made.
    """
    if fig is None:
        fig = plt.figure()

    if ax is None:
        ax = fig.add_axes([0.6, 0, 0.4, 1])

    data = exposure.image.array
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(data)

    ax.imshow(data, origin="lower", vmin=vmin, vmax=vmax)

    nrows = len(donut_stamps)
    for i in range(nrows):
        donut = donut_stamps[i]
        xy = donut.centroid_position

        # plot the cross marking that the donut was used
        ax.scatter(xy[0], xy[1], s=200, marker="+", c="m", lw=4)

        # plot the donut number on the plot
        xtext, ytext = xy[0], xy[1]
        ytext -= 60
        if xtext + 100 > 4096:
            xtext -= 250
        if len(str(i)) > 1:  # move to the left label that is too long
            xtext -= 340
        else:
            xtext -= 260
        ax.text(xtext, ytext, f"{i}", fontsize=17, c="white")
    ax.yaxis.tick_right()
    ax.set_xlabel("x [px]")
    ax.set_ylabel("y [px]")
    ax.yaxis.set_label_position("right")
    info = exposure.getMetadata()
    ax.set_title(f"exposure {info['DAYOBS']}{info['SEQNUM']}")
    return


def plot_donut_stamps(donut_stamps):
    """Plot the donut stamp image cutouts.

    Parameters
    ----------
    donut_stamps: Collection of postage stamps as
        lsst.afw.image.maskedImage.MaskedImage.
        with additional metadata.
    """
    # calculate number of rows given
    # the constraint of the number of
    # columns
    n_donuts = len(donut_stamps)
    ncols = 4
    nrows = n_donuts // ncols
    if nrows * ncols < n_donuts:
        nrows += 1

    fig, axs = plt.subplots(nrows, ncols, figsize=(3 * ncols, 3 * nrows))
    ax = np.ravel(axs)
    for i in range(n_donuts):
        donut = donut_stamps[i]
        ax[i].imshow(donut.stamp_im.image.array, origin="lower")
        ax[i].text(80, 80, f"{i}", fontsize=17, c="white")
    fig.subplots_adjust(hspace=0.35)

    # if there are more axes than donuts,
    # turn off the extra axes
    ncells = nrows * ncols
    if ncells > n_donuts:
        for axis in ax[n_donuts:]:
            axis.axis("off")
    plt.show()
    return


def preview_exposures(
    year_month_day,
    exp_start,
    exp_end,
    dataset_ref_or_type="raw",
    collection="LATISS/raw/all",
    instrument="LATISS",
    detector=0,
):
    """Plot auxTel exposures in a given exposure ID range.

    Parameters:
    ----------
    year_month_day: str
       The year, month, and day of observations as a string,
       eg. '20210908'.
    exp_start: int
       The number of the first exposure in a series to plot, eg. 487.
    exp_end: int
       The number of the final exposure in a series to plot, eg. 490.
    dataset_ref_or_type: str, optional
        Dataset ref or type for butler,
        eg. raw or postISRCCD (default: 'raw').
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
    range(exp_start, exp_end) exist.
    """

    butler = dafButler.Butler("/repo/main/")

    # figure out how many images to plot
    nexp = exp_end - exp_start

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
    for exp in range(exp_start, exp_end):
        exposure = butler.get(
            dataset_ref_or_type,
            dataId={
                "instrument": instrument,
                "detector": detector,
                "exposure": int(f"{year_month_day}00{exp}"),
            },
            collections=[collection],
        )
        data = exposure.image.array
        vmin, vmax = zscale.get_limits(data)
        ax[i].imshow(data, vmin=vmin, vmax=vmax, origin="lower")
        ax[i].set_title(
            f"{year_month_day}, exp {exp},\n focusz={np.round(exposure.getMetadata()['FOCUSZ'],3)}"
        )
        i += 1
    fig.subplots_adjust(hspace=0.35)
    fig.subplots_adjust(wspace=0.3)

    # if there are more axes than exposures,
    # turn off the extra axes
    ncells = nrows * ncol
    if ncells > nexp:
        for axis in ax[nexp:]:
            axis.axis("off")
