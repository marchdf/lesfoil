"""Utilities for scripts."""
import numpy as np
import yaml

airfoil_aoa = 13.3


def parse_ic(fname):
    """Parse the Nalu yaml input file for the initial conditions."""
    with open(fname, "r") as stream:
        try:
            dat = yaml.full_load(stream)
            u0 = float(
                dat["realms"][0]["initial_conditions"][0]["value"]["velocity"][0]
            )
            rho0 = float(
                dat["realms"][0]["material_properties"]["specifications"][0]["value"]
            )
            mu = float(
                dat["realms"][0]["material_properties"]["specifications"][1]["value"]
            )
            turb_model = dat["realms"][0]["solution_options"]["turbulence_model"]
            dt = dat["Time_Integrators"][0]["StandardTimeIntegrator"]["time_step"]

            return u0, rho0, mu, turb_model, dt

        except yaml.YAMLError as exc:
            print(exc)


def ccw_rotation(x, y, angle=airfoil_aoa, rotcenx=0.0, rotceny=0.0, scale=1.0):
    """Rotate counter clockwise."""
    theta = np.radians(angle)
    crdvec = np.array([np.cos(theta), -np.sin(theta)])
    tan_crdvec = np.array([np.sin(theta), np.cos(theta)])
    xp = np.dot(np.asarray([x - rotcenx, y - rotceny]).T, crdvec) / scale + rotcenx
    yp = np.dot(np.asarray([x - rotcenx, y - rotceny]).T, tan_crdvec) / scale + rotceny
    return xp, yp


def ccw_rotation_t00(t00, t01, t11, angle=airfoil_aoa):
    """Return 00 component of rotated tensor."""
    theta = np.radians(angle)
    return (
        t00 * (np.cos(theta) ** 2)
        - t01 * np.sin(2 * theta)
        + t11 * (np.sin(theta) ** 2)
    )


def ccw_rotation_t01(t00, t01, t11, angle=airfoil_aoa):
    """Return 01 component of rotated tensor."""
    theta = np.radians(angle)
    return (
        0.5 * t00 * np.sin(2 * theta)
        + t01 * np.cos(2 * theta)
        - 0.5 * t11 * np.sin(2 * theta)
    )


def ccw_rotation_t11(t00, t01, t11, angle=airfoil_aoa):
    """Return 11 component of rotated tensor."""
    theta = np.radians(angle)
    return (
        t00 * (np.sin(theta) ** 2)
        + t01 * np.sin(2 * theta)
        + t11 * (np.cos(theta) ** 2)
    )


def cord_locations():
    """Return cord locations for slicing."""
    return [0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.825, 0.87, 0.93, 0.99]


def lo_cord_locations():
    """Low cord locations."""
    return [x for x in cord_locations() if x < 0.55]


def hi_cord_locations():
    """High cord locations."""
    return [x for x in cord_locations() if x >= 0.55]


def lo_idx(x, val):
    """Index in x where x <= val."""
    lo = np.argmin(np.fabs(x - val))
    if x[lo] >= val:
        lo = lo - 1
    return lo
