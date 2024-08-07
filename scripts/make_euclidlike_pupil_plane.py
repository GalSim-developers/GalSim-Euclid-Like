"""
This code has been taken from https://github.com/CosmoStat/wf-psf Liaudat et al.

This is an implementation of the function `generate_pupil_obscurations` https://github.com/CosmoStat/wf-psf/blob/87e0c8e9770199cd276f5f0551054cb4902d53bb/src/wf_psf/sims/SimPSFToolkit.py#L233

NOTE from Tobias Liaudat:
"Simple procedure considering only the 2D plane.
No 3D projections wrt the angle of the FoV is done."
"""  # noqa
from argparse import ArgumentParser
import os
import importlib.util
from importlib.resources import files

import numpy as np
from scipy.signal import convolve2d

import galsim


# Telescope parameters
AS_diam = 1200  # Aperture stop diameter [mm]
M1_diam = 395  # Mirror 1 cap stopper diameter [mm]

sp_lenght = 700  # Spider length [mm]
sp_width = 12  # Spider width [mm]

AS_centre = [0, 0]
M1_centre = [0, 51]

sp1_angle = 106.78 - 90  # [degrees]
sp2_angle = 50.11 - 90  # [degrees]
sp3_angle = -10.76 - 90  # [degrees]

sp1_x_pos = 260  # [mm]
sp1_y_pos = 240  # [mm]
sp2_x_pos = -330  # [mm]
sp2_y_pos = 130  # [mm]
sp3_x_pos = 70  # [mm]
sp3_y_pos = -330  # [mm]


def round_even(n):
    return int(2 * np.round(n / 2))


def make_EuclidLike_pupil_plane(N_pix=2048, do_filter=True, N_filter=3):
    # Build pupil plane
    pupil_plane = np.ones((N_pix, N_pix))

    # coordinates of map in [mm]
    W, H = np.meshgrid(
        np.linspace(-AS_diam // 2, AS_diam // 2, N_pix),
        np.linspace(-AS_diam // 2, AS_diam // 2, N_pix),
    )

    # Calculate the Aperture stop and draw it ###
    aperture_stop_mask = np.sqrt(
        (W - AS_centre[0]) ** 2 + (H - AS_centre[1]) ** 2
    ) <= (AS_diam / 2)
    pupil_plane[~aperture_stop_mask] = 0

    # Calculate the M1/M2 obscurations and draw them ###
    M1_mask = np.sqrt((W - M1_centre[0]) ** 2 + (H - M1_centre[1]) ** 2) <= (
        M1_diam / 2
    )
    pupil_plane[M1_mask] = 0

    # Calculate the spiders and draw them
    # Spider 1
    sp1_a = np.tan(sp1_angle * (np.pi / 180))
    sp1_b = sp1_y_pos - sp1_a * sp1_x_pos

    sp1_mask_1 = sp1_a * W + sp1_b - sp_width / 2 * np.sqrt(1 + sp1_a**2) < H
    sp1_mask_2 = sp1_a * W + sp1_b + sp_width / 2 * np.sqrt(1 + sp1_a**2) > H
    sp1_mask = np.logical_and(sp1_mask_1, sp1_mask_2)

    sp1_length_mask = np.sqrt((W - sp1_x_pos) ** 2 + (H - sp1_y_pos) ** 2) <= (
        sp_lenght / 2
    )
    sp1_mask = np.logical_and(sp1_mask, sp1_length_mask)

    # Spider 2
    sp2_a = np.tan(sp2_angle * (np.pi / 180))
    sp2_b = sp2_y_pos - sp2_a * sp2_x_pos

    sp2_mask_1 = sp2_a * W + sp2_b - sp_width / 2 * np.sqrt(1 + sp2_a**2) < H
    sp2_mask_2 = sp2_a * W + sp2_b + sp_width / 2 * np.sqrt(1 + sp2_a**2) > H
    sp2_mask = np.logical_and(sp2_mask_1, sp2_mask_2)

    sp2_length_mask = np.sqrt((W - sp2_x_pos) ** 2 + (H - sp2_y_pos) ** 2) <= (
        sp_lenght / 2
    )
    sp2_mask = np.logical_and(sp2_mask, sp2_length_mask)

    # Spider 3
    sp3_a = np.tan(sp3_angle * (np.pi / 180))
    sp3_b = sp3_y_pos - sp3_a * sp3_x_pos

    sp3_mask_1 = sp3_a * W + sp3_b - sp_width / 2 * np.sqrt(1 + sp3_a**2) < H
    sp3_mask_2 = sp3_a * W + sp3_b + sp_width / 2 * np.sqrt(1 + sp3_a**2) > H
    sp3_mask = np.logical_and(sp3_mask_1, sp3_mask_2)

    sp3_length_mask = np.sqrt((W - sp3_x_pos) ** 2 + (H - sp3_y_pos) ** 2) <= (
        sp_lenght / 2
    )
    sp3_mask = np.logical_and(sp3_mask, sp3_length_mask)

    # Draw the three spider arms
    pupil_plane[sp1_mask] = 0
    pupil_plane[sp2_mask] = 0
    pupil_plane[sp3_mask] = 0

    # Low-pass filter the image
    if do_filter:
        top_hat_filter = np.ones((N_filter, N_filter))

        pupil_plane = convolve2d(
            pupil_plane,
            top_hat_filter,
            boundary="fill",
            mode="same",
            fillvalue=0,
        )
        pupil_plane /= np.sum(top_hat_filter)

    # Get pupil scale
    pupil_scale = AS_diam / N_pix

    return pupil_plane, pupil_scale


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-o",
        "--output",
        dest="output_path",
        help="Path where to save the pupil plane. Default to the euclidlike"
        " data directory if it exist.",
        type=str,
    )
    parser.add_argument(
        "-N",
        "--N_pix",
        dest="Npix",
        help="N pixels used to save the pupil plane. [Default: 2048]",
        type=int,
        default=2048,
    )
    parser.add_argument(
        "-filter",
        "--do_filter",
        dest="DoFilter",
        help="Wether to apply a Top-Hat filter. [Default: True]",
        type=bool,
        default=True,
    )
    parser.add_argument(
        "-Nf",
        "--N_filter",
        dest="Nfilter",
        help="N pixels used for the Top-Hat filter. [Default: 3]",
        type=int,
        default=3,
    )

    args = parser.parse_args()

    # Deal with the output
    if args.output_path is None:
        if importlib.util.find_spec("euclidlike") is None:
            raise ValueError(
                "No output path provided and euclidlike not found."
                " Please provide an output path."
            )
        output_path = str(
            files("euclidlike.data").joinpath("euclid_pupil_plane.fits.gz")
        )
    else:
        output_path = os.path.abspath(args.output_path)

    pupil_plane, pupil_scale = make_EuclidLike_pupil_plane(
        args.Npix,
        args.DoFilter,
        args.Nfilter,
    )

    # Galsim need the pupil scale in m/pixel
    pupil_img = galsim.Image(pupil_plane, scale=pupil_scale / 1000)
    pupil_img.write(output_path)
    print(f"Pupil plane saved at: {output_path}")
