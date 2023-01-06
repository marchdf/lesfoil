#!/usr/bin/env python3

import argparse
import pathlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np
from scipy import interpolate
import utilities


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
        "-f", "--fname", help="Data file to plot", required=True, type=str
    )
    parser.add_argument(
        "-y", "--yname", help="Input file used", required=True, type=str
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
        [RefData("dxi", "fig3a-", x) for x in ["CM1", "CM2", "CM3"]]
        + [RefData("deta", "fig3b-", x) for x in ["CM1", "CM2", "CM3"]]
        + [RefData("dzeta", "fig3c-", x) for x in ["CM1", "CM2", "CM3"]]
        + [RefData("dt", "fig3d-", x) for x in ["CM1", "CM2", "CM3"]]
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
            rd.xdata(), rd.ydata(), **styles[rd.val], label=labels[rd.val],
        )

    # Estimate wall units (utau/nu)
    N_zeta_CM3 = 705
    L_span = 0.0493
    dzeta_physical = L_span / (N_zeta_CM3 - 1)
    rd_zeta = RefData("dzeta", "fig3c-", "CM3")
    wall_units_zeta = rd_zeta.ydata() / dzeta_physical

    rd_eta = RefData("deta", "fig3b-", "CM3")
    eta_interp = interpolate.interp1d(
        rd_zeta.xdata(), rd_zeta.ydata(), bounds_error=False
    )
    deta_physical = rd_eta.ydata() / eta_interp(rd_eta.xdata()) * dzeta_physical
    deta_physical_mean = np.mean(
        deta_physical[np.where((0.3 < rd_eta.xdata()) & (rd_eta.xdata() < 0.7))]
    )
    wall_units_eta = rd_eta.ydata() / deta_physical_mean

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
        rd_zeta.xdata(), wall_units_zeta, lw=2, color=cmap[0], label="First estimate",
    )
    plt.plot(
        rd_eta.xdata(),
        wall_units_eta,
        lw=2,
        color=cmap[1],
        ls="--",
        label="Second estimate",
    )
    plt.plot(
        xnew, wall_units_mean, lw=2, color=cmap[2], ls="--", label="Mean",
    )

    plt.figure("deta_physical")
    plt.plot(
        rd_eta.xdata(), deta_physical, lw=2, color=cmap[0], label=labels[rd_eta.val],
    )
    plt.axhline(
        deta_physical_mean,
        xmin=-100,
        xmax=100,
        lw=2,
        color=cmap[-1],
        ls="--",
        label="Mean",
    )

    # Nalu-Wind data
    u0, rho0, mu, turb_model, dt = utilities.parse_ic(args.yname)
    model = turb_model.upper().replace("_", "-")
    cord = 1.0
    refArea = 0.05
    tau = cord / u0
    dynPres = rho0 * 0.5 * u0 * u0
    re = rho0 * u0 * cord / mu
    aoa = 13.3
    rotcen = 0.0

    cpcf = pd.read_csv(args.fname)
    cpcf["cf"] = cpcf.tauw / dynPres
    cpcf["cfx"] = cpcf.tauwx / dynPres
    cpcf["cfy"] = cpcf.tauwy / dynPres
    cpcf["cp"] = cpcf.pressure / dynPres
    cpcf.sort_values(by=["theta"], inplace=True)
    cord_angle = np.radians(aoa)
    crdvec = np.array([np.cos(cord_angle), -np.sin(cord_angle)])
    tan_crdvec = np.array([np.sin(cord_angle), np.cos(cord_angle)])
    cpcf["xovc"] = (
        np.dot(np.asarray([cpcf.x - rotcen, cpcf.y]).T, crdvec) / cord + rotcen
    )
    cpcf["yovc"] = (
        np.dot(np.asarray([cpcf.x - rotcen, cpcf.y]).T, tan_crdvec) / cord + rotcen
    )
    cpcf["cfxp"] = cpcf.cfx * crdvec[0] + cpcf.cfy * crdvec[1]
    cpcf["cfyp"] = -cpcf.cfx * crdvec[1] + cpcf.cfy * crdvec[0]

    idx_xmax = np.argmax(cpcf.xovc)
    upper = np.where(
        ((cpcf.yovc > 0) & (cpcf.xovc <= 0.5))
        | ((cpcf.yovc > cpcf.yovc.iloc[idx_xmax]) & (cpcf.xovc > 0.5))
    )
    lower = [i for i in range(len(cpcf)) if i not in upper[0].tolist()]

    plt.figure("cp")
    p = plt.plot(cpcf.xovc, -cpcf.cp, lw=2, color=cmap[3], label=f"{model}",)

    plt.figure("airfoil")
    p = plt.plot(cpcf.x, cpcf.y, lw=2, color=cmap[0], label="physical",)
    p = plt.plot(
        cpcf.xovc.iloc[lower],
        cpcf.yovc.iloc[lower],
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

    plt.figure("cf")
    p = plt.plot(
        cpcf.xovc.iloc[upper],
        cpcf.cfx.iloc[upper],
        lw=2,
        color=cmap[3],
        label=f"{model}",
    )

    # Save the plots
    fname = "plots.pdf"
    with PdfPages(fname) as pdf:

        plt.figure("dxi")
        ax = plt.gca()
        plt.xlabel(r"$x/c$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$\Delta \xi^+$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        plt.xlim([0, 1])
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("deta")
        ax = plt.gca()
        plt.xlabel(r"$x/c$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$\Delta \eta^+$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        plt.xlim([0, 1])
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("dzeta")
        ax = plt.gca()
        plt.xlabel(r"$x/c$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$\Delta \zeta^+$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        plt.xlim([0, 1])
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("dt")
        ax = plt.gca()
        plt.xlabel(r"$x/c$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$\Delta t^+$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        plt.xlim([0, 1])
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("wall_units")
        ax = plt.gca()
        plt.xlabel(r"$x/c$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$u_\tau / \nu$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        plt.xlim([0, 1])
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("deta_physical")
        ax = plt.gca()
        plt.xlabel(r"$x/c$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$\Delta \eta$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        plt.xlim([0, 1])
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("cp")
        ax = plt.gca()
        plt.xlabel(r"$x/c$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$-C_p$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        plt.xlim([0, 1])
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("cf")
        ax = plt.gca()
        plt.xlabel(r"$x/c$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$C_f$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        plt.xlim([0, 1])
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("u-a")
        ax = plt.gca()
        plt.xlabel(r"$\langle u \rangle / u_\infty$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$y_\eta / c$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        plt.xlim([0, 7])
        plt.ylim([0, None])
        # remove duplicate labels
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        legend = ax.legend(by_label.values(), by_label.keys(), loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("u-b")
        ax = plt.gca()
        plt.xlabel(r"$\langle u \rangle / u_\infty$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$y_\eta / c$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        plt.xlim([0, 7])
        plt.ylim([0, None])
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        legend = ax.legend(by_label.values(), by_label.keys(), loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("urms-a")
        ax = plt.gca()
        plt.xlabel(r"$\langle u' \rangle / u_\infty$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$y_\eta / c$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        plt.xlim([-0.1, 1.5])
        plt.ylim([0, None])
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        legend = ax.legend(by_label.values(), by_label.keys(), loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("urms-b")
        ax = plt.gca()
        plt.xlabel(r"$\langle u' \rangle / u_\infty$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$y_\eta / c$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        plt.xlim([-0.1, 1.5])
        plt.ylim([0, None])
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        legend = ax.legend(by_label.values(), by_label.keys(), loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("vrms-a")
        ax = plt.gca()
        plt.xlabel(r"$\langle v' \rangle / u_\infty$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$y_\eta / c$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        plt.xlim([-0.1, 1.5])
        plt.ylim([0, None])
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        legend = ax.legend(by_label.values(), by_label.keys(), loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("vrms-b")
        ax = plt.gca()
        plt.xlabel(r"$\langle v' \rangle / u_\infty$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$y_\eta / c$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        plt.xlim([-0.1, 1.5])
        plt.ylim([0, None])
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        legend = ax.legend(by_label.values(), by_label.keys(), loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("uvrms-a")
        ax = plt.gca()
        plt.xlabel(r"$\langle u'v' \rangle / u_\infty$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$y_\eta / c$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        plt.xlim([-0.01, 0.07])
        plt.ylim([0, None])
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        legend = ax.legend(by_label.values(), by_label.keys(), loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("uvrms-b")
        ax = plt.gca()
        plt.xlabel(r"$\langle u'v' \rangle / u_\infty$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$y_\eta / c$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        plt.xlim([-0.01, 0.07])
        plt.ylim([0, None])
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
