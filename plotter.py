#!/usr/bin/env python3

import argparse
import pathlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd


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


if __name__ == "__main__":

    # Parse arguments
    parser = argparse.ArgumentParser(description="A simple plot tool")
    args = parser.parse_args()

    # Reference data
    refdir = pathlib.Path("refdata")

    
    # rdf = pd.read_csv(os.path.join(refdir, "data.txt"))
    # print(rdf)
    # tmdir = os.path.join(refdir, "tomboulides2000")
    # ngdir = os.path.join(refdir, "nagata2020")

    dct = {"dxi":{"pfx":"fig3a-", "vals": ["CM1", "CM2", "CM3"]}}

    for key,v in dct.items():

        for i, val in enumerate(v["vals"]):
            fn = v["pfx"] + val + ".csv"
            df = pd.read_csv(refdir / fn, header=None, names=["x", key])
            print(df)
            plt.figure(key)
            plt.plot(df.x, df[key], lw=2, color=cmap[i], label=f"""Asada \& Kawai (2018) {val}""")
    
    # # Plot reference values
    # grouped = rdf.groupby("ref")
    # for k, (name, group) in enumerate(grouped):
    #     for coef in ["cd", "cl"]:
    #         plt.figure(f"{coef}-re")
    #         p = plt.plot(
    #             group.Re,
    #             group[coef],
    #             lw=0,
    #             color=cmap[-1],
    #             linestyle="None",
    #             marker=markertype[k],
    #             ms=10,
    #             label=name,
    #         )
    #         plt.figure(f"{coef}-spectra")
    #         for cnt, st in enumerate(group.st):
    #             if not np.isnan(st):
    #                 Re = group.iloc[cnt].Re
    #                 p = plt.axvline(
    #                     st,
    #                     ymin=0,
    #                     ymax=30,
    #                     lw=2,
    #                     color=cmap[-1],
    #                     ls="--",
    #                     label=f"{name}, Re = {Re}",
    #                 )

    # # Nagata specific
    # plt.figure("cp")
    # df = pd.read_csv(os.path.join(ngdir, "cp.csv"))
    # p = plt.plot(
    #     df.theta,
    #     savgol_filter(df.cp, 21, 4),
    #     lw=2,
    #     color=cmap[-1],
    #     label="Nagata (2020), Re=300, Ma=0.3",
    # )
    # # p = plt.plot(df.theta, df.cp, lw=2, color=cmap[-1])
    # plt.figure("cf")
    # df = pd.read_csv(os.path.join(ngdir, "cf.csv"))
    # p = plt.plot(
    #     df.theta,
    #     savgol_filter(df.cf, 21, 4),
    #     lw=2,
    #     color=cmap[-1],
    #     label="Nagata (2020), Re=300, Ma=0.3",
    # )
    # # p = plt.plot(df.theta, df.cf, lw=2, color=cmap[-1])

    # # Tomboulides specific
    # plt.figure("cp")
    # df = pd.read_csv(os.path.join(tmdir, "cp.csv"))
    # p = plt.plot(
    #     df.theta,
    #     savgol_filter(df.cp, 21, 4),
    #     lw=2,
    #     color=cmap[-1],
    #     label="Tomboulides (2000), Re=100",
    # )
    # # p = plt.plot(df.theta, df.cp, lw=2, color=cmap[-1])

    # # AMR-Wind/Nalu-Wind data
    # lst = []
    # for i, fdir in enumerate(args.fdir):

    #     yname = os.path.join(fdir, "sphere-nwind.yaml")
    #     u0, rho0, mu, turb_model, dt = utilities.parse_ic(yname)
    #     dynPres = rho0 * 0.5 * u0 * u0
    #     diameter = 1.0
    #     refArea = np.pi * diameter * diameter / 4.0
    #     tau = diameter / u0
    #     Re = round(u0 / mu)
    #     model = turb_model.upper().replace("_", "-")
    #     ename = os.path.join(fdir, "exwsim.yaml")
    #     is_hybrid = os.path.exists(ename)
    #     if is_hybrid:
    #         sim = "AMR-Wind + Nalu-Wind"
    #         rdir = "overset"
    #     else:
    #         sim = "Nalu-Wind"
    #         rdir = "results"
    #     lbl = f"{sim}, Re={Re}"
    #     legend_elements += [
    #         Line2D([0], [0], lw=2, color=cmap[i], label=lbl),
    #     ]

    #     all_files = glob.glob(os.path.join(fdir, "forces*.dat"))
    #     df = pd.concat(
    #         (pd.read_csv(f, delim_whitespace=True) for f in all_files)
    #     ).sort_values(by="Time", ignore_index=True)
    #     df["t"] = df.Time / tau
    #     df = df[df["t"] > args.tmin]
    #     df["px"] = df.Fpx / (dynPres * refArea)
    #     df["vx"] = df.Fvx / (dynPres * refArea)
    #     df["cd"] = (df.Fpx + df.Fvx) / (dynPres * refArea)
    #     df["cy"] = (df.Fpy + df.Fvy) / (dynPres * refArea)
    #     df["cz"] = (df.Fpz + df.Fvz) / (dynPres * refArea)
    #     alpha = utilities.find_symmetry_plane(df)
    #     df["cl"] = -(np.cos(alpha) * df.cy - np.sin(alpha) * df.cz)
    #     if df.cl.mean() > 0.0:
    #         df.cl *= -1
    #     df["cs"] = np.sin(alpha) * df.cy + np.cos(alpha) * df.cz

    #     means = df.mean()
    #     stats = {}
    #     stats["Re"] = Re
    #     for coef in ["cd", "cl", "cs"]:
    #         stats[coef] = means[coef]
    #         stats[coef + "_err"] = df[coef].max() - means[coef]
    #         stats[coef + "_amp"] = df[coef].max() - df[coef].min()

    #         plt.figure(coef)
    #         p = plt.plot(df["t"], df[coef], lw=2, color=cmap[i], label=lbl)
    #         p[0].set_dashes(dashseq[i])

    #         fhat, yhat = utilities.getFFT(dt, df[coef], normalize=True, window=False)
    #         yhat = np.abs(yhat[1:])
    #         fhat = fhat[1:]
    #         stats[coef + "_st"] = fhat[np.argmax(yhat)]
    #         plt.figure(coef + "-spectra")
    #         p = plt.plot(fhat * tau, yhat, lw=2, color=cmap[i], label=lbl)

    #     wdf = pd.read_csv(os.path.join(fdir, rdir, "wall.dat"))
    #     wdf.sort_values(by=["theta"], inplace=True)
    #     wdf["cp"] = wdf.pressure / dynPres
    #     wdf["cf"] = wdf.tauwx / dynPres

    #     bins = np.linspace(0, 180, 100)
    #     group = wdf.groupby(pd.cut(wdf.theta, bins))
    #     plot_centers = (bins[:-1] + bins[1:]) / 2
    #     cp_values = group.cp.mean()
    #     cf_values = group.cf.mean()

    #     plt.figure("cp")
    #     # p = plt.plot(wdf.theta, wdf.cp, lw=0, color=cmap[-1], linestyle="None", marker="s", ms=1)
    #     p = plt.plot(plot_centers, cp_values, lw=2, color=cmap[i], label=lbl)
    #     p[0].set_dashes(dashseq[i])

    #     plt.figure("cf")
    #     # p = plt.plot(wdf.theta, wdf.cf, lw=0, color=cmap[-1], linestyle="None", marker="s", ms=1)
    #     p = plt.plot(plot_centers, cf_values, lw=2, color=cmap[i], label=lbl)
    #     p[0].set_dashes(dashseq[i])

    #     print(f"""Stats for {fdir}""")
    #     print(f"""Symmetry plane is at alpha = {np.degrees(alpha)}""")
    #     print(stats)
    #     lst.append(stats)

    # df = pd.DataFrame(lst)
    # plt.figure("cd-re")
    # p = plt.plot(
    #     df.Re,
    #     df.cd,
    #     lw=0,
    #     color=cmap[0],
    #     linestyle="None",
    #     marker="s",
    #     ms=10,
    #     mew=2,
    #     mfc="none",
    #     label="AMR-Wind + Nalu-Wind",
    # )
    # plt.figure("cl-re")
    # p = plt.plot(
    #     df.Re,
    #     df.cl,
    #     lw=0,
    #     color=cmap[0],
    #     linestyle="None",
    #     marker="s",
    #     ms=10,
    #     mfc="none",
    #     mew=2,
    #     label="AMR-Wind + Nalu-Wind",
    # )

    # Save the plots
    fname = "plots.pdf"
    with PdfPages(fname) as pdf:

        plt.figure("dxi")
        ax = plt.gca()
        plt.xlabel(r"$x/c$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$\Delta \xi^+$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        # plt.xlim([0, 30])
        #plt.ylim([0.4, 1.1])
        # plt.ylim([0.6, 0.7])
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        # plt.figure("cl")
        # ax = plt.gca()
        # plt.xlabel(r"$t$", fontsize=22, fontweight="bold")
        # plt.ylabel(r"$c_L$", fontsize=22, fontweight="bold")
        # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        # # plt.xlim([0, 30])
        # # plt.ylim([-0.1, 0])
        # legend = ax.legend(loc="best")
        # plt.tight_layout()
        # pdf.savefig(dpi=300)

        # plt.figure("cs")
        # ax = plt.gca()
        # plt.xlabel(r"$t$", fontsize=22, fontweight="bold")
        # plt.ylabel(r"$c_S$", fontsize=22, fontweight="bold")
        # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        # # plt.xlim([0, 30])
        # plt.ylim([-0.1, 0.1])
        # legend = ax.legend(loc="best")
        # plt.tight_layout()
        # pdf.savefig(dpi=300)

        # plt.figure("cd-spectra")
        # ax = plt.gca()
        # plt.xlabel(r"$St$", fontsize=22, fontweight="bold")
        # plt.ylabel(r"$\hat{C}_D$", fontsize=22, fontweight="bold")
        # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        # plt.xlim([0, 1])
        # plt.ylim([0, 2e-3])
        # legend = ax.legend(loc="best")
        # plt.tight_layout()
        # pdf.savefig(dpi=300)

        # plt.figure("cl-spectra")
        # ax = plt.gca()
        # plt.xlabel(r"$St$", fontsize=22, fontweight="bold")
        # plt.ylabel(r"$\hat{C}_L$", fontsize=22, fontweight="bold")
        # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        # plt.xlim([0, 1])
        # plt.ylim([0, 1e-2])
        # legend = ax.legend(loc="best")
        # plt.tight_layout()
        # pdf.savefig(dpi=300)

        # plt.figure("cp")
        # ax = plt.gca()
        # plt.xlabel(r"$\theta$", fontsize=22, fontweight="bold")
        # plt.ylabel(r"$c_P$", fontsize=22, fontweight="bold")
        # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        # # plt.xlim([0, 30])
        # # plt.ylim([0.6, 0.7])
        # legend = ax.legend(loc="best")
        # plt.tight_layout()
        # pdf.savefig(dpi=300)

        # plt.figure("cf")
        # ax = plt.gca()
        # plt.xlabel(r"$\theta$", fontsize=22, fontweight="bold")
        # plt.ylabel(r"$c_f$", fontsize=22, fontweight="bold")
        # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        # # plt.xlim([0, 30])
        # # plt.ylim([0.6, 0.7])
        # legend = ax.legend(loc="best")
        # plt.tight_layout()
        # pdf.savefig(dpi=300)

        # plt.figure("cd-re")
        # ax = plt.gca()
        # plt.xlabel(r"$Re$", fontsize=22, fontweight="bold")
        # plt.ylabel(r"$c_D$", fontsize=22, fontweight="bold")
        # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        # # plt.xlim([0, 30])
        # # plt.ylim([0.4, 1.1])
        # # plt.ylim([0.6, 0.7])
        # legend = ax.legend(loc="best")
        # plt.tight_layout()
        # pdf.savefig(dpi=300)

        # plt.figure("cl-re")
        # ax = plt.gca()
        # plt.xlabel(r"$Re$", fontsize=22, fontweight="bold")
        # plt.ylabel(r"$c_L$", fontsize=22, fontweight="bold")
        # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        # # plt.xlim([0, 30])
        # # plt.ylim([0.4, 1.1])
        # # plt.ylim([0.6, 0.7])
        # legend = ax.legend(loc="best")
        # plt.tight_layout()
        # pdf.savefig(dpi=300)
