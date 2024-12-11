"""
This script processes Euclid observation sequence XML file and convert it into
a pandas DataFrame.
It handle the conversion from ecliptic coordinates (from EC) to equatorial.
The dates of the observations are derived from Figure 4.4 of MOCDC_v4.2
(https://euclid.roe.ac.uk/dmsf/files/20821/view  EC only).

NOTE: The true observation date is shifted by 6 years to match the current
Roman/Rubin simulations. That allows us to include transients.

NOTE 2: This observation sequence has been made from the `rsd2024a` 
https://euclid.roe.ac.uk/dmsf/files/20839/view (EC only).

To run the script run the command:
    python make_obseq.py <path_to_xml_file>
"""
import os
import sys
from tqdm import tqdm
import xml.etree.ElementTree as ET
from argparse import ArgumentParser
import importlib.util
from importlib.resources import files

import pandas as pd
import numpy as np

from euclidlike.instrument_params import (
    long_exptime,
    short_exptime_vis,
    nisp_exptime_total,
    nisp_exptime_eff,
)

from astropy.time import Time
from astropy.coordinates import SkyCoord, GeocentricTrueEcliptic
import astropy.units as u


def s2days(s):
    return s / 86400.


def convert_ecliptic_to_equatorial(elon, elat, epa):
    """ Transform coordinate frame.

    Parameters
    ----------
    elon : np.ndarray
        longitude coordinate (deg)
    elat : np.ndarray
        latitude coordinate (deg)
    epa : np.ndarray
        position angle (deg)

    Returns
    -------
    lon, lat, pa : np.ndarray
    """

    skycoord = SkyCoord(
                elon*u.deg,
                elat*u.deg,
                frame=GeocentricTrueEcliptic,
    )

    out = skycoord.transform_to("icrs")

    lon = out.ra.deg
    lat = out.dec.deg

    pole_target = SkyCoord(
        0*u.deg,
        90*u.deg,
        frame="icrs"
    )

    pole_0 = pole_target.transform_to(GeocentricTrueEcliptic)

    pa = (epa - skycoord.position_angle(pole_0).deg) % 360

    return lon, lat, pa


shift_2_years = 365 * 2 + 1
shift_6_years = 365 * 6 + 1

VIS_exp_long = long_exptime
VIS_exp_short = short_exptime_vis
NISP_exp_int = nisp_exptime_total
NISP_exp_eff = nisp_exptime_eff
NISP_exp_wait = NISP_exp_int - NISP_exp_eff

# Figure 4.4 of MOCDC_v4.2
# https://euclid.roe.ac.uk/dmsf/files/20821/view (EC only)
dither_seq = {
    0: {
        "VIS": [
            {"wait": 8},
            {"VIS_LONG": VIS_exp_long},
            {"wait": 40 + NISP_exp_int + 22 + NISP_exp_int + 22},
            {"VIS_SHORT": VIS_exp_short},
        ],
        "NISP": [
            {"wait": 574 + 40},
            {"NISP_J": NISP_exp_eff},
            {"wait": 22 + NISP_exp_wait},
            {"NISP_H": NISP_exp_eff},
            {"wait": 22 + NISP_exp_wait},
            {"NISP_Y": NISP_exp_eff},
        ],
    },
    1: {
        "VIS": [
            {"wait": 8},
            {"VIS_LONG": VIS_exp_long},
            {"wait": 40 + NISP_exp_int + 22 + NISP_exp_int + 22},
            {"VIS_SHORT": VIS_exp_short},
        ],
        "NISP": [
            {"wait": 574 + 40},
            {"NISP_J": NISP_exp_eff},
            {"wait": 22 + NISP_exp_wait},
            {"NISP_H": NISP_exp_eff},
            {"wait": 22 + NISP_exp_wait},
            {"NISP_Y": NISP_exp_eff},
        ],
    },
    2: {
        "VIS": [
            {"wait": 8},
            {"VIS_LONG": VIS_exp_long},
            {"wait": 40 + NISP_exp_int + 22 + NISP_exp_int + 22},
            {"VIS_SHORT": 0},
        ],
        "NISP": [
            {"wait": 574 + 40},
            {"NISP_J": NISP_exp_eff},
            {"wait": 22 + NISP_exp_wait},
            {"NISP_H": NISP_exp_eff},
            {"wait": 22 + NISP_exp_wait},
            {"NISP_Y": NISP_exp_eff},
        ],
    },
    3: {
        "VIS": [
            {"wait": 8},
            {"VIS_LONG": VIS_exp_long},
            {"wait": 40 + NISP_exp_int + 22 + NISP_exp_int + 22},
            {"VIS_SHORT": 0},
        ],
        "NISP": [
            {"wait": 574 + 40},
            {"NISP_J": NISP_exp_eff},
            {"wait": 22 + NISP_exp_wait},
            {"NISP_H": NISP_exp_eff},
            {"wait": 22 + NISP_exp_wait},
            {"NISP_Y": NISP_exp_eff},
        ],
    },
}

def make_obseq(input_path, output_path):
    date = []
    date_euclid = []
    obs_id = []
    pointing_id = []
    dither_id = []
    center_ra = []
    center_dec = []
    pa = []
    center_lon = []
    center_lat = []
    pa_ecl = []
    saa = []

    tree = ET.parse(input_path)
    root = tree.getroot()
    obs_req = list(root.iter("ObservationRequest"))
    for obs in tqdm(obs_req, total=len(obs_req)):
        survey_name =  obs.find('SurveyId').text
        if survey_name != "WIDE_SURVEY":
            continue

        obs_id_ = int(obs.get('id'))

        for pointing in obs.iter('PointingRequest'):
            # Convert from mjd2000 to mjd
            mjd = int(pointing.find('Mjd2000').text)/100_000 + 51544.0
            saa_ = float(pointing.find("SAA").text)
            attitude = pointing.find('Attitude')
            lon_ = float(attitude.find('Longitude').text)
            lat_ = float(attitude.find('Latitude').text)
            pa_ = float(attitude.find('PositionAngle').text)
            pointing_id_ = int(pointing.attrib['id'])
            dither_id_ = int(pointing.attrib['ditherId'])

            t = Time(mjd, format="mjd")
            new_t = t.mjd + shift_6_years

            # Ecliptic to ra, dec
            ra_, dec_, pa_eq_ = convert_ecliptic_to_equatorial(lon_, lat_, pa_)

            date.append(new_t)
            date_euclid.append(t.mjd)
            obs_id.append(obs_id_)
            pointing_id.append(pointing_id_)
            dither_id.append(dither_id_)
            center_ra.append(ra_)
            center_dec.append(dec_)
            pa.append(pa_eq_)
            center_lon.append(lon_)
            center_lat.append(lat_)
            pa_ecl.append(pa_)
            saa.append(saa_)

    col_dtype = [
        ("date", np.float64),
        ("exptime", np.float64),
        ("ra", np.float64),
        ("dec", np.float64),
        ("pa", np.float64),
        ("lon", np.float64),
        ("lat", np.float64),
        ("pa_ecliptic", np.float64),
        ("saa", np.float64),
        ("filter", "<U6"),
        ("date_euclid", np.float64),
        ("obs_id", np.int64),
        ("pointing_id", np.int64),
        ("patch_id", np.int64),
        ("dither_id", np.int16),
    ]
    obs_kind = {
        "VIS_LONG": 0,
        "VIS_SHORT": 1,
        "NISP_J": 2,
        "NISP_H": 3,
        "NISP_Y": 4,
    }

    data = np.zeros((len(date), len(obs_kind)), dtype=col_dtype)

    for i in tqdm(range(len(date)), total=len(date)):
        dither = dither_id[i]
        for instrum in ["VIS", "NISP"]:
            for j, obs_kind_ in enumerate(dither_seq[dither][instrum]):
                obs_kind_name = list(obs_kind_.keys())[0]
                if j == 0:
                    d = date[i]
                    d_eucl = date_euclid[i]
                if obs_kind_name == "wait":
                    d += s2days(obs_kind_["wait"])
                    d_eucl += s2days(obs_kind_["wait"])
                    continue
                if "VIS" in obs_kind_name:
                    filter_name = "VIS"
                else:
                    filter_name = obs_kind_name
                obs_kind_ind = obs_kind[obs_kind_name]
                data["date"][i, obs_kind_ind] = d
                data["exptime"][i, obs_kind_ind] = obs_kind_[obs_kind_name]
                data["ra"][i, obs_kind_ind] = center_ra[i]
                data["dec"][i, obs_kind_ind] = center_dec[i]
                data["pa"][i, obs_kind_ind] = pa[i]
                data["lon"][i, obs_kind_ind] = center_lon[i]
                data["lat"][i, obs_kind_ind] = center_lat[i]
                data["pa_ecliptic"][i, obs_kind_ind] = pa_ecl[i]
                data["saa"][i, obs_kind_ind] = saa[i]
                data["filter"][i, obs_kind_ind] = filter_name
                data["date_euclid"][i, obs_kind_ind] = d_eucl
                data["obs_id"][i, obs_kind_ind] = obs_id[i]
                data["pointing_id"][i, obs_kind_ind] = pointing_id[i]
                data["dither_id"][i, obs_kind_ind] = dither

    df = pd.DataFrame(
        data.reshape(len(date)*len(obs_kind)),
        index=pd.MultiIndex.from_product(
            [
                range(len(date)),
                obs_kind,
            ],
            names=["visit", "obs_kind"],
        ),
    )
    df.to_pickle(output_path)


def parse_args(command_args):
    parser = ArgumentParser()
    parser.add_argument(
        "input_xml",
        help="Input obseq file from Euclid (.xml format).",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_path",
        help="Path where to save the obseq file. Default to the euclidlike_imsim"
        " data directory if it exist.",
        type=str,
    )

    args = parser.parse_args(command_args)
    return args


def main(command_args):

    args = parse_args(command_args)

    # Deal with the output
    if args.output_path is None:
        if importlib.util.find_spec("euclidlike_imsim") is None:
            raise ValueError(
                "No output path provided and euclidlike_imsim not found."
                " Please provide an output path."
            )
        output_path = str(
            files("euclidlike_imsim.data").joinpath("euclid_obseq.pkl")
        )
    else:
        output_path = os.path.abspath(args.output_path)

    make_obseq(args.input_xml, output_path)

    print(f"obseq saved at: {output_path}")


def run_main():
    main(sys.argv[1:])
