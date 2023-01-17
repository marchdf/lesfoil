#!/usr/bin/env python3

import argparse
import glob
import os
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import utilities as ut
from matplotlib import ticker
from matplotlib.backends.backend_pdf import PdfPages
from scipy import interpolate

plt.rc("text", usetex=True)
plt.rcParams.update({"figure.max_open_warning": 0})
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


class RefData:
    """Reference data class"""

    def __init__(self, key, pfx, val, sort=None, sfx=None):
        """Initialize."""
        self.key = key
        self.pfx = pfx
        self.val = val
        self.sort = sort
        self.sfx = sfx
        if self.sfx:
            self.fname = self.pfx + self.val + self.sfx + ".csv"
        else:
            self.fname = self.pfx + self.val + ".csv"
        self.df = pd.read_csv(refdir / self.fname, header=None, names=["x", key])
        self.df.sort_values(by="x", inplace=True)

        if self.sort == "polar":
            center = [0.1, 0.0]
            self.df["angles"] = np.arctan2(
                self.df.x - center[0], self.df[key] - center[1]
            )
            self.df.sort_values(by="angles", inplace=True)
        elif self.sort == "y":
            self.df.sort_values(by=key, inplace=True)

    def xdata(self):
        """Return xdata."""
        return self.df["x"].to_numpy()

    def ydata(self):
        """Return ydata."""
        return self.df[self.key].to_numpy()


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description="A simple plot tool")
    parser.add_argument(
        "-f", "--fdir", help="Folder with data files", required=True, type=str
    )
    args = parser.parse_args()

    # Reference data
    refdir = pathlib.Path("refdata")
    styles = {
        "CM1": {"lw": 2, "color": cmap[0]},
        "CM2": {"lw": 2, "color": cmap[1]},
        "CM3": {"lw": 2, "color": cmap[2]},
        "exp": {
            "lw": 0,
            "color": cmap[-1],
            "linestyle": "None",
            "marker": markertype[0],
            "mfc": "None",
            "ms": 6,
        },
        "lesfoil": {
            "lw": 0,
            "color": cmap[-1],
            "linestyle": "None",
            "marker": markertype[1],
            "mfc": "None",
            "ms": 6,
        },
    }
    labels = {
        "CM1": "Asada \& Kawai (2018) CM1",
        "CM2": "Asada \& Kawai (2018) CM2",
        "CM3": "Asada \& Kawai (2018) CM3",
        "exp": "Gleyzes (1978) Exp.",
        "lesfoil": "Mary \& Sagaut (2002) LES",
    }
    rdata = (
        [RefData("dxip", "fig3a-", x) for x in ["CM1", "CM2", "CM3"]]
        + [RefData("detap", "fig3b-", x) for x in ["CM1", "CM2", "CM3"]]
        + [RefData("dzetap", "fig3c-", x) for x in ["CM1", "CM2", "CM3"]]
        + [RefData("dtp", "fig3d-", x) for x in ["CM1", "CM2", "CM3"]]
        + [
            RefData("cp", "fig4a-", x, sort="polar")
            for x in ["CM1", "CM2", "CM3", "exp", "lesfoil"]
        ]
        + [RefData("cf", "fig4b-", x) for x in ["CM1", "CM2", "CM3", "exp", "lesfoil"]]
        + [RefData("u-a", "fig5a-", "CM3", sfx=f"-{x}", sort="y") for x in range(5)]
        + [RefData("u-a", "fig5a-", "exp", sfx=f"-{x}", sort="y") for x in range(3, 5)]
        + [RefData("u-a", "fig5a-", "lesfoil", sfx=f"-{x}", sort="y") for x in range(5)]
        + [RefData("u-b", "fig5b-", "CM3", sfx=f"-{x}", sort="y") for x in range(5)]
        + [RefData("u-b", "fig5b-", "exp", sfx=f"-{x}", sort="y") for x in range(5)]
        + [RefData("u-b", "fig5b-", "lesfoil", sfx=f"-{x}", sort="y") for x in range(5)]
        + [RefData("urms-a", "fig6a-", "CM3", sfx=f"-{x}", sort="y") for x in range(5)]
        + [
            RefData("urms-a", "fig6a-", "exp", sfx=f"-{x}", sort="y")
            for x in range(3, 5)
        ]
        + [
            RefData("urms-a", "fig6a-", "lesfoil", sfx=f"-{x}", sort="y")
            for x in range(5)
        ]
        + [RefData("urms-b", "fig6b-", "CM3", sfx=f"-{x}", sort="y") for x in range(5)]
        + [RefData("urms-b", "fig6b-", "exp", sfx=f"-{x}", sort="y") for x in range(5)]
        + [
            RefData("urms-b", "fig6b-", "lesfoil", sfx=f"-{x}", sort="y")
            for x in range(5)
        ]
        + [RefData("vrms-a", "fig7a-", "CM3", sfx=f"-{x}", sort="y") for x in range(5)]
        + [
            RefData("vrms-a", "fig7a-", "exp", sfx=f"-{x}", sort="y")
            for x in range(3, 5)
        ]
        + [RefData("vrms-b", "fig7b-", "CM3", sfx=f"-{x}", sort="y") for x in range(5)]
        + [RefData("vrms-b", "fig7b-", "exp", sfx=f"-{x}", sort="y") for x in range(5)]
        + [RefData("uvrms-a", "fig8a-", "CM3", sfx=f"-{x}", sort="y") for x in range(5)]
        + [
            RefData("uvrms-a", "fig8a-", "exp", sfx=f"-{x}", sort="y")
            for x in range(3, 5)
        ]
        + [
            RefData("uvrms-a", "fig8a-", "lesfoil", sfx=f"-{x}", sort="y")
            for x in range(5)
        ]
        + [RefData("uvrms-b", "fig8b-", "CM3", sfx=f"-{x}", sort="y") for x in range(5)]
        + [RefData("uvrms-b", "fig8b-", "exp", sfx=f"-{x}", sort="y") for x in range(5)]
        + [
            RefData("uvrms-b", "fig8b-", "lesfoil", sfx=f"-{x}", sort="y")
            for x in range(5)
        ]
    )

    for rd in rdata:
        plt.figure(rd.key)
        plt.plot(
            rd.xdata(),
            rd.ydata(),
            **styles[rd.val],
            label=labels[rd.val],
        )

    # Estimate wall units (utau/nu)
    N_zeta_CM3 = 705
    L_span = 0.0493
    dzeta = L_span / (N_zeta_CM3 - 1)
    rd_zeta = RefData("dzetap", "fig3c-", "CM3")
    wall_units_zeta = rd_zeta.ydata() / dzeta

    rd_eta = RefData("detap", "fig3b-", "CM3")
    eta_interp = interpolate.interp1d(
        rd_zeta.xdata(), rd_zeta.ydata(), bounds_error=False
    )
    deta = rd_eta.ydata() / eta_interp(rd_eta.xdata()) * dzeta
    deta_mean = np.mean(deta[np.where((0.3 < rd_eta.xdata()) & (rd_eta.xdata() < 0.7))])
    wall_units_eta = rd_eta.ydata() / deta_mean

    wall_units_zeta_interp = interpolate.interp1d(
        rd_zeta.xdata(), wall_units_zeta, bounds_error=False
    )
    wall_units_eta_interp = interpolate.interp1d(
        rd_eta.xdata(), wall_units_eta, bounds_error=False
    )
    xnew = np.linspace(0, 1, 1000)
    wall_units_mean = 0.5 * (wall_units_zeta_interp(xnew) + wall_units_eta_interp(xnew))
    wall_units_mean_interp = interpolate.interp1d(
        xnew, wall_units_mean, bounds_error=False
    )

    plt.figure("wall_units")
    plt.plot(
        rd_zeta.xdata(),
        wall_units_zeta,
        lw=2,
        color=cmap[0],
        label=labels[rd_zeta.val] + " estimate",
    )
    plt.plot(
        rd_eta.xdata(),
        wall_units_eta,
        lw=2,
        color=cmap[1],
        ls="--",
        label=labels[rd_zeta.val] + " estimate",
    )
    plt.plot(
        xnew,
        wall_units_mean,
        lw=2,
        color=cmap[2],
        ls="--",
        label=labels[rd_zeta.val] + " mean",
    )

    rd_dxip = RefData("dxip", "fig3a-", "CM3")
    plt.figure("dxi")
    plt.plot(
        rd_dxip.xdata(),
        rd_dxip.ydata() / wall_units_mean_interp(rd_dxip.xdata()),
        lw=2,
        color=cmap[0],
        label=labels[rd_dxip.val] + " estimate",
    )
    plt.figure("deta")
    plt.plot(
        rd_eta.xdata(),
        deta,
        lw=2,
        color=cmap[0],
        label=labels[rd_eta.val] + " estimate",
    )
    plt.axhline(
        deta_mean,
        xmin=-100,
        xmax=100,
        lw=2,
        color=cmap[-1],
        ls="--",
        label="Mean",
    )
    plt.figure("dzeta")
    plt.axhline(
        dzeta,
        xmin=-100,
        xmax=100,
        lw=2,
        color=cmap[0],
        label=labels[rd_zeta.val],
    )
    rd_dtp = RefData("dtp", "fig3d-", "CM3")
    plt.figure("dt")
    plt.plot(
        rd_dtp.xdata(),
        rd_dtp.ydata() / wall_units_mean_interp(rd_dtp.xdata()),
        lw=2,
        color=cmap[0],
        label=labels[rd_dxip.val] + " estimate",
    )

    # Nalu-Wind data
    yname = glob.glob(os.path.join(args.fdir, "*.yaml"))[0]
    fname = os.path.join(args.fdir, "wing.dat")
    pname = os.path.join(args.fdir, "profiles.dat")
    u0, rho0, mu, turb_model, dt = ut.parse_ic(yname)
    nu = mu / rho0
    model = turb_model.upper().replace("_", "-")
    cord = 1.0
    refArea = 0.05
    tau = cord / u0
    dynPres = rho0 * 0.5 * u0 * u0
    re = rho0 * u0 * cord / mu
    deta = 0.000011813977015662547  # from the PW mesh

    # wing data
    cpcf = pd.read_csv(fname)
    cpcf["cf"] = cpcf.tauw / dynPres
    cpcf["cfx"] = cpcf.tauwx / dynPres
    cpcf["cfy"] = cpcf.tauwy / dynPres
    cpcf["cp"] = cpcf.pressure / dynPres
    cpcf.sort_values(by=["theta"], inplace=True)
    cpcf["xovc"], cpcf["yovc"] = ut.ccw_rotation(cpcf.x, cpcf.y)
    cpcf["ref_wall_units"] = wall_units_mean_interp(cpcf.xovc)
    cpcf["dx"] = np.diff(cpcf.x, append=cpcf.x.iloc[0])
    cpcf["dy"] = np.diff(cpcf.y, append=cpcf.y.iloc[0])
    cpcf["dxi"] = np.sqrt(cpcf.dx**2 + cpcf.dy**2)
    cpcf["deta"] = deta
    cpcf["dzeta"] = cpcf.dz
    cpcf["dt"] = dt
    cpcf["dxip"] = cpcf.dxi * cpcf.ref_wall_units
    cpcf["detap"] = cpcf.deta * cpcf.ref_wall_units
    cpcf["dzetap"] = cpcf.dzeta * cpcf.ref_wall_units
    cpcf["dtp"] = cpcf.dt * cpcf.ref_wall_units
    cpcf["utau"] = np.sqrt(cpcf.tauw / rho0)
    cpcf["wall_units"] = cpcf.utau / nu
    tgt = np.vstack(
        (np.diff(cpcf.x, append=cpcf.x.iloc[0]), np.diff(cpcf.y, append=cpcf.y.iloc[0]))
    )
    tgt /= np.linalg.norm(tgt, axis=0)
    cpcf["cftgt"] = cpcf.cfx * tgt[0, :] + cpcf.cfy * tgt[1, :]
    cpcf["cfnml"] = -cpcf.cfx * tgt[1, :] + cpcf.cfy * tgt[0, :]

    idx_xmax = np.argmax(cpcf.xovc)
    upper = np.where(
        ((cpcf.yovc > 0) & (cpcf.xovc <= 0.5))
        | ((cpcf.yovc > cpcf.yovc.iloc[idx_xmax]) & (cpcf.xovc > 0.5))
    )
    lower = [i for i in range(len(cpcf)) if i not in upper[0].tolist()]

    for val in [
        "dxip",
        "detap",
        "dzetap",
        "dtp",
        "dxi",
        "deta",
        "dzeta",
        "dt",
        "wall_units",
    ]:
        plt.figure(val)
        plt.plot(
            cpcf.xovc.iloc[upper],
            cpcf[val].iloc[upper],
            lw=2,
            color=cmap[3],
            label=f"{model}",
        )

    plt.figure("cp")
    p = plt.plot(
        cpcf.xovc,
        -cpcf.cp,
        lw=2,
        color=cmap[3],
        label=f"{model}",
    )

    plt.figure("cf")
    p = plt.plot(
        cpcf.xovc.iloc[upper],
        cpcf.cftgt.iloc[upper],
        lw=2,
        color=cmap[3],
        label=f"{model}",
    )

    plt.figure("airfoil")
    p = plt.plot(
        cpcf.x,
        cpcf.y,
        lw=2,
        color=cmap[0],
        label="physical",
    )
    p = plt.plot(
        cpcf.iloc[lower].sort_values(by=["xovc"]).xovc,
        cpcf.iloc[lower].sort_values(by=["xovc"]).yovc,
        lw=2,
        color=cmap[1],
        label="rotated, lower",
    )
    p = plt.plot(
        cpcf.xovc.iloc[upper],
        cpcf.yovc.iloc[upper],
        lw=2,
        color=cmap[2],
        label="rotated, upper",
    )

    # profiles
    ndf = pd.read_csv(pname)
    ndf["urms"] = ndf.upup - ndf.tau_xx
    ndf["vrms"] = ndf.vpvp - ndf.tau_yy
    ndf["uvrms"] = ndf.upvp - ndf.tau_xy
    grouped = ndf.groupby(["xloc"])
    fields = {
        "u": 1.4,
        "urms": 0.3,
        "vrms": 0.3,
        "uvrms": 0.014,
        "tke": 0.03,
        "sdr": 0.0,
        "tvisc": 0.0004,
        "kratio": 1.0,
        "alpha": 1.0,
    }
    for k, (xloc, group) in enumerate(grouped):
        sfx = "a" if xloc in ut.lo_cord_locations() else "b"
        cnt = (
            ut.lo_cord_locations().index(xloc)
            if xloc in ut.lo_cord_locations()
            else ut.hi_cord_locations().index(xloc)
        )

        for field, offset in fields.items():
            plt.figure(f"{field}-{sfx}")
            p = plt.plot(
                group[field] + cnt * offset,
                group.eta,
                lw=2,
                color=cmap[3],
                label=f"{model}",
            )
            p[0].set_dashes(dashseq[0])

    # Save the plots
    fname = "plots.pdf"
    with PdfPages(fname) as pdf:
        plots = {
            "dxip": {
                "ylabel": r"$\Delta \xi^+$",
            },
            "detap": {
                "ylabel": r"$\Delta \eta^+$",
            },
            "dzetap": {
                "ylabel": r"$\Delta \zeta^+$",
            },
            "dtp": {
                "ylabel": r"$\Delta t^+$",
            },
            "wall_units": {
                "ylabel": r"$u_\tau / \nu$",
            },
            "dxi": {
                "ylabel": r"$\Delta \xi$",
            },
            "deta": {
                "ylabel": r"$\Delta \eta$",
            },
            "dzeta": {
                "ylabel": r"$\Delta \zeta$",
            },
            "dt": {
                "ylabel": r"$\Delta t$",
            },
            "cp": {
                "ylabel": r"$-C_p$",
            },
            "cf": {
                "ylabel": r"$C_f$",
            },
        }
        for name, plot in plots.items():
            plt.figure(name)
            ax = plt.gca()
            plt.xlabel(r"$x/c$", fontsize=22, fontweight="bold")
            plt.ylabel(plot["ylabel"], fontsize=22, fontweight="bold")
            plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
            plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
            plt.xlim([0, 1])
            legend = ax.legend(loc="best")
            plt.tight_layout()
            pdf.savefig(dpi=300)

        sci_format = ticker.ScalarFormatter(useMathText=True)
        sci_format.set_scientific(True)
        sci_format.set_powerlimits((-1, 1))

        plots = {
            "u": {"xlabel": r"$\langle u \rangle / u_\infty$", "xlim": [0, 7]},
            "urms": {
                "xlabel": r"$\langle u' \rangle_{rms} / u_\infty$",
                "xlim": [-0.1, 1.5],
            },
            "vrms": {
                "xlabel": r"$\langle v' \rangle_{rms} / u_\infty$",
                "xlim": [-0.1, 1.5],
            },
            "uvrms": {
                "xlabel": r"$\langle u'v' \rangle / u_\infty^2$",
                "xlim": [-0.01, 0.07],
            },
            "tke": {"xlabel": r"$\langle k \rangle$", "xlim": [-0.01, 0.15]},
            "sdr": {
                "xlabel": r"$\langle \omega \rangle$",
                "xlim": [1, 1e5],
                "xlog": True,
            },
            "tvisc": {
                "xlabel": r"$\langle \nu_t \rangle$",
                "xlim": [-0.0001, 0.002],
                "x_sci_format": sci_format,
            },
            "kratio": {
                "xlabel": r"$\langle k_r \rangle$",
                "xlim": [-0.01, 5],
            },
            "alpha": {
                "xlabel": r"$\langle \alpha \rangle$",
                "xlim": [-0.01, 5],
            },
        }
        for name, plot in plots.items():
            for sub in ["a", "b"]:
                plt.figure(f"{name}-{sub}")
                ax = plt.gca()
                plt.xlabel(plot["xlabel"], fontsize=22, fontweight="bold")
                plt.ylabel(r"$y_\eta / c$", fontsize=22, fontweight="bold")
                plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
                plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
                if "xlog" in plot:
                    ax.set_xscale("log")
                plt.xlim(plot["xlim"])
                plt.ylim([0, 0.012] if sub == "a" else [0, 0.09])
                if "x_sci_format" in plot:
                    ax.xaxis.set_major_formatter(plot["x_sci_format"])
                # remove duplicate labels
                handles, labels = ax.get_legend_handles_labels()
                by_label = dict(zip(labels, handles))
                legend = ax.legend(by_label.values(), by_label.keys(), loc="best")
                plt.tight_layout()
                pdf.savefig(dpi=300)

        plt.figure("airfoil")
        ax = plt.gca()
        plt.xlabel(r"$x/c$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$y / c$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        legend = ax.legend(loc="best")
        ax.axis("equal")
        plt.tight_layout()
        pdf.savefig(dpi=300)
