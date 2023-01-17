"""Post process."""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import stk
import utilities as ut
from matplotlib.backends.backend_pdf import PdfPages
from mpi4py import MPI
from scipy import interpolate
from scipy.interpolate import griddata

plt.rc("text", usetex=True)
cmap_med = [
    "#F15A60",
    "#7AC36A",
    "#5A9BD4",
    "#FAA75B",
    "#9E67AB",
    "#CE7058",
    "#D77FB4",
    "#737373",
]
cmap = [
    "#EE2E2F",
    "#008C48",
    "#185AA9",
    "#F47D23",
    "#662C91",
    "#A21D21",
    "#B43894",
    "#010202",
]
dashseq = [
    (None, None),
    [10, 5],
    [10, 4, 3, 4],
    [3, 3],
    [10, 4, 3, 4, 3, 4],
    [3, 3],
    [3, 3],
]
markertype = ["s", "d", "o", "p", "h"]


def p0_printer(par):
    iproc = par.rank

    def printer(*args, **kwargs):
        if iproc == 0:
            print(*args, **kwargs)

    return printer


def subset_fields(data, xloc, yloc, m_airfoil, radius=0.1):
    # get the coords without the aoa rotation
    xp, yp = ut.ccw_rotation(data[:, 0], data[:, 1], angle=ut.airfoil_aoa())

    # take the data above the cord vector of the aifoil and in a radius around (xloc, yloc)
    return data[
        (yp >= m_airfoil * xp) & (((xp - xloc) ** 2 + (yp - yloc) ** 2) < radius**2),
        :,
    ]


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description="A simple post-processing tool")
    parser.add_argument(
        "-m",
        "--mfile",
        help="Root name of files to postprocess",
        required=True,
        type=str,
    )
    parser.add_argument("--auto_decomp", help="Auto-decomposition", action="store_true")
    parser.add_argument(
        "-v",
        "--vel_name",
        help="Name of the velocity field",
        default="velocity",
        type=str,
    )
    parser.add_argument(
        "--navg", help="Number of times to average", default=40, type=int
    )
    parser.add_argument(
        "--flowthrough", help="Flowthrough time (L/u)", default=1.0, type=float
    )
    parser.add_argument(
        "--factor",
        help="Factor of flowthrough time between time steps used in average",
        type=float,
        default=1.2,
    )
    args = parser.parse_args()

    fdir = os.path.dirname(args.mfile)
    mname = os.path.splitext(os.path.basename(args.mfile))[0]

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    par = stk.Parallel.initialize()
    printer = p0_printer(par)

    odir = os.path.join(os.getcwd(), f"output-{mname}")
    if rank == 0:
        if not os.path.exists(odir):
            os.makedirs(odir)

    mesh = stk.StkMesh(par)
    printer("Reading meta data for mesh: ", args.mfile)
    mesh.read_mesh_meta_data(args.mfile, auto_decomp=args.auto_decomp)
    printer("Done reading meta data")

    printer("Loading bulk data for mesh: ", args.mfile)
    mesh.populate_bulk_data()
    printer("Done reading bulk data")

    num_time_steps = mesh.stkio.num_time_steps
    max_time = mesh.stkio.max_time
    tsteps = np.array(mesh.stkio.time_steps)
    printer(f"""Num. time steps = {num_time_steps}""")
    printer(f"""Max. time step  = {max_time}""")
    printer(f"""All steps in file: {tsteps}""")

    # Figure out the times over which to average
    if args.factor > 0:
        tmp_tavg = np.sort(
            tsteps[-1] - args.flowthrough * args.factor * np.arange(args.navg)
        )
        tmp_tavg = tmp_tavg[tmp_tavg >= tsteps[0]]
        dist = np.abs(np.array(tsteps)[:, np.newaxis] - tmp_tavg)
        idx = dist.argmin(axis=0)
    else:
        idx = np.arange(len(tsteps) - args.navg, len(tsteps))
    tavg = tsteps[idx[0] :]
    tavg_instantaneous = tsteps[idx[0] :]
    # printer("Averaging the following steps:")
    # printer(tavg)

    # Extract time and spanwise average tau_wall on wall
    wingdata = None
    for tstep in tavg_instantaneous:
        ftime, missing = mesh.stkio.read_defined_input_fields(tstep)
        printer(f"Loading tau_wall fields for time: {ftime}")

        coords = mesh.meta.coordinate_field
        wing = mesh.meta.get_part("airfoil")
        sel = wing & mesh.meta.locally_owned_part
        tauw = mesh.meta.get_field("tau_wall")
        tauwv = mesh.meta.get_field("tau_wall_vector")
        pressure = mesh.meta.get_field("pressure")
        names = ["x", "y", "z", "tauw", "tauwx", "tauwy", "tauwz", "pressure"]
        nnodes = sum(bkt.size for bkt in mesh.iter_buckets(sel, stk.StkRank.NODE_RANK))

        cnt = 0
        data = np.zeros((nnodes, len(names)))
        for bkt in mesh.iter_buckets(sel, stk.StkRank.NODE_RANK):
            xyz = coords.bkt_view(bkt)
            tw = tauw.bkt_view(bkt)
            twv = tauwv.bkt_view(bkt)
            pres = pressure.bkt_view(bkt)
            data[cnt : cnt + bkt.size, :] = np.hstack(
                (xyz, tw.reshape(-1, 1), twv, pres.reshape(-1, 1))
            )
            cnt += bkt.size

        if wingdata is None:
            wingdata = np.zeros(data.shape)
        wingdata += data / len(tavg_instantaneous)

    lst = comm.gather(wingdata, root=0)
    comm.Barrier()
    if rank == 0:
        df = pd.DataFrame(np.vstack(lst), columns=names)
        dz = np.diff(np.unique(df.z)).min()
        df = df.groupby("x", as_index=False).mean().sort_values(by=["x"])
        center = [0.1, 0.0]
        df["theta"] = np.arctan2(df.x - center[0], df.y - center[1])
        df["dz"] = dz
        wingname = os.path.join(odir, "wing.dat")
        df.to_csv(wingname, index=False)

    # Extract velocity profiles
    is_ams = not mesh.meta.get_field("average_velocity").is_null
    vel_name = "velocity"
    dudx_name = "dudx"
    field_names = [
        "u",
        "v",
        "w",
        "tke",
        "sdr",
        "tvisc",
        "kratio",
        "alpha",
        "tau_xx",
        "tau_xy",
        "tau_yy",
    ]
    fld_data = None
    for tstep in tavg:
        ftime, missing = mesh.stkio.read_defined_input_fields(tstep)
        printer(f"""Loading fields for time: {ftime}""")

        interior = mesh.meta.get_part("flow-hex")
        sel = interior & mesh.meta.locally_owned_part
        coords = mesh.meta.coordinate_field
        turbulent_ke = mesh.meta.get_field("turbulent_ke")
        specific_dissipation_rate = mesh.meta.get_field("specific_dissipation_rate")
        tvisc = mesh.meta.get_field("turbulent_viscosity")
        fields = [
            mesh.meta.get_field(vel_name),
            turbulent_ke,
            specific_dissipation_rate,
            tvisc,
        ]
        dveldx = mesh.meta.get_field(dudx_name)
        k_ratio = mesh.meta.get_field("k_ratio")
        names = ["x", "y", "z"] + field_names
        nnodes = sum(bkt.size for bkt in mesh.iter_buckets(sel, stk.StkRank.NODE_RANK))

        cnt = 0
        data = np.zeros((nnodes, len(names)))
        for bkt in mesh.iter_buckets(sel, stk.StkRank.NODE_RANK):
            arr = coords.bkt_view(bkt)
            for fld in fields:
                vals = fld.bkt_view(bkt)
                if len(vals.shape) == 1:  # its a scalar
                    vals = vals.reshape(-1, 1)
                arr = np.hstack((arr, vals))

            # tauSGRS_ij = coeffSGRS *(avgdudx[:, i * 3 + j] + avgdudx[:, j * 3 + i]) + 2/3 rho k delta_ij
            dudx = dveldx.bkt_view(bkt)
            nut = tvisc.bkt_view(bkt)
            tke = turbulent_ke.bkt_view(bkt)
            if is_ams:
                alpha = k_ratio.bkt_view(bkt) ** 1.7
                krat = k_ratio.bkt_view(bkt)
            else:
                alpha = np.ones(nut.shape)
                krat = np.ones(nut.shape)
            rho = 1.0
            coeffSGRS = alpha * (2.0 - alpha) * nut / rho
            diag_tke = (-2.0 / 3.0 * rho * tke * krat).reshape(-1, 1)
            tausgrs_xx = (coeffSGRS * (dudx[:, 0] + dudx[:, 0])).reshape(
                -1, 1
            ) + diag_tke
            tausgrs_xy = (coeffSGRS * (dudx[:, 1] + dudx[:, 3])).reshape(-1, 1)
            tausgrs_yy = (coeffSGRS * (dudx[:, 4] + dudx[:, 4])).reshape(
                -1, 1
            ) + diag_tke
            arr = np.hstack((arr, krat.reshape(-1, 1)))
            arr = np.hstack((arr, alpha.reshape(-1, 1)))
            arr = np.hstack((arr, tausgrs_xx))
            arr = np.hstack((arr, tausgrs_xy))
            arr = np.hstack((arr, tausgrs_yy))
            data[cnt : cnt + bkt.size, :] = arr
            cnt += bkt.size

        if fld_data is None:
            fld_data = np.zeros(data.shape)
        fld_data += data / len(tavg)

    # Load airfoil shape
    upper = pd.read_csv(
        "./meshes/lesfoil_upper.dat",
        header=None,
        names=["x", "y", "z"],
        skiprows=1,
        delim_whitespace=True,
    )
    upper.sort_values(by=["x"])
    upper_y_interp = interpolate.interp1d(upper.x, upper.y, bounds_error=False)
    m_airfoil = (upper.y.iloc[-1] - upper.y.iloc[0]) / (
        upper.x.iloc[-1] - upper.x.iloc[0]
    )

    comm.Barrier()
    if rank == 0:
        plt.figure("airfoil")
        p = plt.plot(
            upper.x,
            upper.y,
            lw=2,
            color="red",
            label="upper",
        )

    # Subset the fields
    ninterp = 200
    deta = 0.095
    planes = []
    for xloc in ut.cord_locations():
        yloc = upper_y_interp(xloc)

        sub = subset_fields(data, xloc, yloc, m_airfoil)

        lst = comm.gather(sub, root=0)
        comm.Barrier()
        if rank == 0:
            # normal to upper part of the airfoil
            lo_idx = ut.lo_idx(upper.x.to_numpy(), xloc)
            tgt = (upper.iloc[lo_idx + 1] - upper.iloc[lo_idx]).to_numpy()
            tgt /= np.linalg.norm(tgt)
            nml = [-tgt[1], tgt[0], 0]

            # equation of normal
            m = nml[1] / nml[0]
            p = yloc - m * xloc
            dxnml = np.sqrt(deta**2 / (1 + m**2))
            xnml = (
                np.linspace(xloc - dxnml, xloc, ninterp)
                if m < 0
                else np.linspace(xloc, xloc + dxnml, ninterp)
            )
            ynml = m * xnml + p

            plt.figure("airfoil")
            p = plt.plot(
                xloc,
                yloc,
                lw=0,
                color="black",
                linestyle="None",
                marker="o",
                mfc="None",
                ms=6,
            )
            p = plt.plot(xnml, ynml, lw=2, color="green", label="normal")

            df = (
                pd.DataFrame(np.vstack(lst), columns=names)
                .groupby(["x", "y"], as_index=False)
                .mean()
                .sort_values(by=["x", "y"])
            )

            # rotate the data to remove the aoa rotation
            df["xa"], df["ya"] = ut.ccw_rotation(df.x, df.y)
            df["ua"], df["va"] = ut.ccw_rotation(df.u, df.v)
            df["tau_xxa"] = ut.ccw_rotation_t00(df.tau_xx, df.tau_xy, df.tau_yy)
            df["tau_xya"] = ut.ccw_rotation_t01(df.tau_xx, df.tau_xy, df.tau_yy)
            df["tau_yya"] = ut.ccw_rotation_t11(df.tau_xx, df.tau_xy, df.tau_yy)

            # rotate data so that the tangent is horizontal and the normal is vertical
            angle = np.degrees(np.arctan2(tgt[1], tgt[0]))
            df["x"], df["y"] = ut.ccw_rotation(df.xa - xloc, df.ya - yloc, angle=-angle)
            df["u"], df["v"] = ut.ccw_rotation(df.ua, df.va, angle=-angle)
            df["tau_xx"] = ut.ccw_rotation_t00(
                df.tau_xx, df.tau_xy, df.tau_yy, angle=-angle
            )
            df["tau_xy"] = ut.ccw_rotation_t01(
                df.tau_xx, df.tau_xy, df.tau_yy, angle=-angle
            )
            df["tau_yy"] = ut.ccw_rotation_t11(
                df.tau_xx, df.tau_xy, df.tau_yy, angle=-angle
            )
            xp = (upper.x - xloc) * tgt[0] + (upper.y - yloc) * tgt[1]
            yp = -(upper.x - xloc) * tgt[1] + (upper.y - yloc) * tgt[0]
            xi = np.array([0])
            yi = np.logspace(-5, np.log10(deta), ninterp)
            plt.figure(f"airfoil-{xloc}")
            plt.tripcolor(df.x, df.y, df.u, shading="gouraud")
            p = plt.plot(xi * np.ones(yi.shape), yi, lw=2, color="red")
            p = plt.plot(xp, yp, lw=2, color="green", label="Rotated airfoil")
            p = plt.plot(
                0.0,
                0.0,
                lw=0,
                color="black",
                linestyle="None",
                marker="o",
                mfc="None",
                ms=6,
            )

            means = {}
            for fld in field_names:
                means[fld] = griddata(
                    (df.x, df.y),
                    df[fld],
                    (xi[None, :], yi[:, None]),
                    method="cubic",
                    fill_value=0,
                ).flatten()
            means["xloc"] = xloc * np.ones(yi.shape)
            means["yloc"] = yloc * np.ones(yi.shape)
            means["eta"] = yi
            planes.append(pd.DataFrame(means))

    printer("Extract fluctuating velocities")
    comm.Barrier()
    if rank == 0:
        for plane in planes:
            plane["upup"] = np.zeros(plane.u.shape)
            plane["vpvp"] = np.zeros(plane.u.shape)
            plane["upvp"] = np.zeros(plane.u.shape)

    # for tstep in tavg_instantaneous:
    #     ftime, missing = mesh.stkio.read_defined_input_fields(tstep)
    #     printer(f"""Loading {vel_name} fields for time: {ftime}""")

    #     interior = mesh.meta.get_part("flow-hex")
    #     sel = interior & mesh.meta.locally_owned_part
    #     coords = mesh.meta.coordinate_field
    #     velocity = mesh.meta.get_field("velocity")
    #     names = ["x", "y", "z", "u", "v", "w"]
    #     nnodes = sum(bkt.size for bkt in mesh.iter_buckets(sel, stk.StkRank.NODE_RANK))

    #     cnt = 0
    #     data = np.zeros((nnodes, len(names)))
    #     for bkt in mesh.iter_buckets(sel, stk.StkRank.NODE_RANK):
    #         xyz = coords.bkt_view(bkt)
    #         vel = velocity.bkt_view(bkt)
    #         data[cnt : cnt + bkt.size, :] = np.hstack((xyz, vel))
    #         cnt += bkt.size

    #     for k, xloc in enumerate(ut.cord_locations()):
    #         sub = subset_fields(data, xloc, yloc, m_airfoil)

    #         lst = comm.gather(sub, root=0)
    #         comm.Barrier()
    #         if rank == 0:
    #             xi = np.array([0])
    #             yi = np.logspace(-5, np.log10(deta), ninterp)
    #             df = pd.DataFrame(np.vstack(lst), columns=names)

    #             # rotate the data to remove the aoa rotation
    #             df["xa"], df["ya"] = ut.ccw_rotation(df.x, df.y, angle=ut.airfoil_aoa())
    #             df["ua"], df["va"] = ut.ccw_rotation(df.u, df.v, angle=ut.airfoil_aoa())

    #             # rotate data so that the tangent is horizontal and the normal is vertical
    #             df["x"] = (df.xa - xloc) * tgt[0] + (df.ya - yloc) * tgt[1]
    #             df["y"] = -(df.xa - xloc) * tgt[1] + (df.ya - yloc) * tgt[0]
    #             df["u"] = df.ua * tgt[0] + df.va * tgt[1]
    #             df["v"] = -df.ua * tgt[1] + df.va * tgt[0]
    #             df.drop(columns=["xa", "ya", "ua", "va"])

    #             grouped = df.groupby("z")
    #             navg = len(tavg_instantaneous) * grouped.ngroups
    #             for name, group in grouped:
    #                 up = (
    #                     griddata(
    #                         (group.x, group.y),
    #                         group.u,
    #                         (xi[None, :], yi[:, None]),
    #                         method="cubic",
    #                         fill_value=0,
    #                     ).flatten()
    #                     - planes[k].u
    #                 )
    #                 vp = (
    #                     griddata(
    #                         (group.x, group.y),
    #                         group.v,
    #                         (xi[None, :], yi[:, None]),
    #                         method="cubic",
    #                         fill_value=0,
    #                     ).flatten()
    #                     - planes[k].v
    #                 )

    #                 planes[k].upup += np.sqrt(up * up) / navg
    #                 planes[k].vpvp += np.sqrt(vp * vp) / navg
    #                 planes[k].upvp += up * vp / navg

    comm.Barrier()
    if rank == 0:
        df = pd.concat(planes)
        df.to_csv(os.path.join(odir, "profiles.dat"), index=False)

        fname = "pp.pdf"
        with PdfPages(fname) as pdf:
            plt.figure("airfoil")
            ax = plt.gca()
            plt.xlabel(r"$x/c$", fontsize=22, fontweight="bold")
            plt.ylabel(r"$y / c$", fontsize=22, fontweight="bold")
            plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
            plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
            # remove duplicate labels
            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            legend = ax.legend(by_label.values(), by_label.keys(), loc="best")
            ax.axis("equal")
            plt.tight_layout()
            pdf.savefig(dpi=300)

            for xloc in ut.cord_locations():
                plt.figure(f"airfoil-{xloc}")
                ax = plt.gca()
                plt.xlabel(r"$(x-x_c)' / c$", fontsize=22, fontweight="bold")
                plt.ylabel(r"$(y-y_c)' / c$", fontsize=22, fontweight="bold")
                plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
                plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
                legend = ax.legend(loc="best")
                ax.axis("equal")
                plt.tight_layout()
                pdf.savefig(dpi=300)
