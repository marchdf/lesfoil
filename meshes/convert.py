"""Convert UFR geometry to something I can use in PW."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


def resample(df):
    """Resample."""
    lst = []
    npt = 10
    for i in range(len(df.x) - 1):
        lst.append(
            np.linspace(df.x.iloc[i], df.x.iloc[i + 1], npt, endpoint=False)
        )
    lst.append(np.array(df.x.iloc[-1]))
    xp = np.hstack(lst)
    f = interp1d(df.x, df.y)
    yp = f(xp)

    dff = pd.DataFrame({"x": xp, "y": yp})
    dff["z"] = 0.0
    # dff["xc"] = dff.x - dff.x.mean()
    # dff["yc"] = dff.y - dff.y.mean()
    # dff["r"] = np.sqrt(dff.xc ** 2 + dff.yc ** 2)
    # dff["theta"] = (np.arctan2(dff.xc, dff.yc) + np.pi * 0.5) * 180 / np.pi
    # dff.sort_values(by=["theta"], inplace=True)
    return dff


def write(df, fname):
    """Write the dataframe."""
    df.to_csv(
        fname, sep=" ", index=False, header=False, columns=["x", "y", "z"]
    )
    with open(fname, "r+") as f:
        content = f.read()
        f.seek(0, 0)
        f.write(str(len(df)) + "\n" + content)


if __name__ == "__main__":
    gname = "UFR2-05_geometry.dat"
    fname = "lesfoil.dat"
    uname = "lesfoil_upper.dat"
    lname = "lesfoil_lower.dat"
    data = []
    with open(gname, "r") as f:
        for line in f:
            if line.startswith("%"):
                continue
            data += [float(x) for x in line.split()]

    df = pd.DataFrame({"x": data[0:513], "y": data[513:]})
    df["z"] = 0.0

    # remove data points not on airfoil
    df = df[df.x <= 1.0]

    # lower and upper
    x0, y0 = 0.0, 0.0
    x1, y1 = 1.0, -0.011846
    m = (y1 - y0) / (x1 - x0)
    p = y0 - m * x0
    ydivider = m * df.x.values + p
    df_upper = df[df.y >= ydivider]
    df_lower = df[df.y <= ydivider]

    # # resample
    # dff_upper = resample(df_upper)
    # dff_lower = resample(df_lower)

    # write file
    df.sort_values(by=["x"], inplace=True)
    write(df, fname)
    write(df_upper, uname)
    write(df_lower, lname)

    plt.figure()
    ax = plt.gca()
    plt.plot(df_upper.x, df_upper.y, "o", ms=1)
    # plt.plot(dff_upper.x, dff_upper.y, "s", ms=0.1)
    plt.plot(df_lower.x, df_lower.y, "o", ms=1)
    # plt.plot(dff_lower.x, dff_lower.y, "s", ms=0.1)
    ax.axis("equal")
    # ax.set(xlim=(0, 0.05), ylim=(-0.05, 0.05))
    plt.savefig("lesfoil.png", dpi=300)
