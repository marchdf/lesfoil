"""Utilities for scripts."""
import yaml


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
