"""Utilities for scripts."""
import yaml
import numpy as np


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


def airfoil_aoa():
    return 13.3


def ccw_rotation(x, y, angle=13.3, rotcenx=0.0, rotceny=0.0, scale=1.0):
    """Rotate counter clockwise."""
    theta = np.radians(angle)
    crdvec = np.array([np.cos(theta), -np.sin(theta)])
    tan_crdvec = np.array([np.sin(theta), np.cos(theta)])
    xp = np.dot(np.asarray([x - rotcenx, y - rotceny]).T, crdvec) / scale + rotcenx
    yp = np.dot(np.asarray([x - rotcenx, y - rotceny]).T, tan_crdvec) / scale + rotceny
    return xp, yp


def cord_locations():
    return [0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.825, 0.87, 0.93, 0.99]


def lo_idx(x, val):
    lo = np.argmin(np.fabs(x - val))
    if x[lo] >= val:
        lo = lo - 1
    return lo
