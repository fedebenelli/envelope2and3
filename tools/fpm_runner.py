import os
from glob import glob

import matplotlib.pyplot as plt

import numpy as np

import pandas as pd


OUTDIR = "env23out/"
LOG = False


def list_outfiles():
    return glob("env23out/envelout*")


def crossed_point(x_values, y_values, point):
    x_point, y_point = point

    dx = x_point - x_values[0]
    dy = y_point - y_values[0]

    for i, (x, y) in enumerate(zip(x_values, y_values)):
        dx_2 = x_point - x
        dy_2 = y_point - y

        if (np.sign(dx) != np.sign(dx_2)) and (np.sign(dy) != np.sign(dy_2)):
            crossed = True
            return i, crossed

        dx = dx_2
        dy = dy_2

    return i, False


def clean_dew(df_dew, dsps):

    t = df_dew["T"].values
    p = df_dew["P"].values

    crossing_index = []

    for dsp in dsps:
        index, crossed = crossed_point(t, p, dsp)

        if crossed:
            crossing_index.append(index)

    crossing_index = sorted(crossing_index)

    if len(crossing_index) == 2:
        i, j = crossing_index
        df_dew[i:j] = np.nan

    elif len(crossing_index) == 1:
        if LOG: print(crossing_index)
        df_dew[crossing_index[0]:] = np.nan

    return df_dew


def clean_bub(df_bub, dsps):
    [t1, p1], [t2, p2] = dsps

    # Clean everything at the left of the lower dsp

    stable = []
    is_stable = False
    crossed_first = False
    crossed_second = False

    for t, p in zip(df_bub["T"], df_bub["P"]):
        if (t - t1) > 0 and (p - p1) > 0 and not crossed_first:
            crossed_first = True

        elif (t - t2) > 0 and (p - p2) < 0 and not crossed_second:
            crossed_second = True

        if crossed_first and not crossed_second:
            is_stable = True
        else:
            is_stable = False

        stable.append(is_stable)

    stable = np.array(stable)
    df_bub[~stable] = np.nan

    return df_bub


def get_stable(envelopes, dsp, hpll):

    if not dsp:
        envelopes["2phase"][1] = np.nan * envelopes["2phase"][1]
    else:
        # DSPs at the start of 3phase lines
        dsp_1 = envelopes["3phase"][0]["T"][0], envelopes["3phase"][0]["P"][0]
        dsp_2 = envelopes["3phase"][2]["T"][0], envelopes["3phase"][2]["P"][0]

        if hpll:
            msk = envelopes["2phase"][2]["T"] > dsp_1[0]
            envelopes["2phase"][2][msk] = np.nan

        envelopes["2phase"][1] = clean_bub(envelopes["2phase"][1], [dsp_1, dsp_2])

        envelopes["2phase"][0] = clean_dew(envelopes["2phase"][0], [dsp_1, dsp_2])

    return envelopes


def set_composition(infile, index=None, value=None, z=None):
    with open(infile) as f:
        lines = f.readlines()

    if index and value:
        z = np.array(lines[1].split(), dtype="float64")
        z[index] = value
        z = z / z.sum()
    elif z is not None:
        pass
    else:
        z = np.array(lines[1].split(), dtype="float64")

    n = len(z)

    lines[0] = str(n) + "\n"
    lines[1] = " ".join(z.astype(str).tolist()) + "\n"
    text = "".join(lines)

    with open("new_file", "w") as w:
        w.write(text)

    return z


def run(infile, three_phase="yes", **kwargs):
    for file in glob(OUTDIR + "*"):
        os.remove(file)

    extra_args = " ".join([f"--{kwarg} {kwargs[kwarg]}" for kwarg in kwargs])
    run_string = (
        f"fpm run {extra_args}"
        f" -- {infile} {three_phase}"
        f" > env23log 2> env23err; wait"
    )
    if LOG: print(run_string)
    os.system(run_string)


def read_envel(file):
    return pd.read_csv(file, delim_whitespace=True)


def read_2phase(pri):
    df_dew = read_envel(OUTDIR + "envelout2-DEW")
    df_bub = read_envel(OUTDIR + "envelout2-LTBUB")

    if pri:
        df_hpll = read_envel(OUTDIR + "envelout2-HPLL")
    else:
        df_hpll = None

    return df_dew, df_bub, df_hpll


def show(df, index, prop, **kwargs):
    if index is None:
        plt.plot(df.index, df[prop], **kwargs)
    else:
        plt.plot(df[index], df[prop], **kwargs)
    return


def get_case():
    outfiles = list_outfiles()
    long_string = " ".join(outfiles)
    if LOG: print(long_string)

    three_phase = True if "envelout3" in long_string else False
    high_pressure_liquid = True if "HPLL" in long_string else False
    dsp_presence = False if "NOCROSS" in long_string else True

    return high_pressure_liquid, dsp_presence, three_phase


def read_nodsp(prim=False):
    df_bub = read_envel(OUTDIR + "envelout2-LTBUB")
    df_dew = read_envel(OUTDIR + "envelout2-DEW")

    df_31 = read_envel(OUTDIR + "envelout3-NOCROSS1")
    df_32 = read_envel(OUTDIR + "envelout3-NOCROSS2")

    return df_bub, df_dew, df_31, df_32


def read_dsp(prim=False):
    df_bub = read_envel(OUTDIR + "envelout2-LTBUB")
    df_dew = read_envel(OUTDIR + "envelout2-DEW")

    df_c11 = read_envel(OUTDIR + "envelout3-CROSS11")
    df_c12 = read_envel(OUTDIR + "envelout3-CROSS12")

    df_c21 = read_envel(OUTDIR + "envelout3-CROSS21")
    df_c22 = read_envel(OUTDIR + "envelout3-CROSS22")

    if prim:
        df_hpll = read_envel(OUTDIR + "envelout2-HPLL")
    else:
        df_hpll = None

    return df_bub, df_dew, df_c11, df_c12, df_c21, df_c22, df_hpll


def read_envels():
    prim, dsp, three_phase = get_case()

    if dsp:
        df_bub, df_dew, df_c11, df_c12, df_c21, df_c22, df_hpll = read_dsp(prim)
        return {
            "prima_case": prim,
            "has_dsp": dsp,
            "2phase": [df_dew, df_bub, df_hpll],
            "3phase": [df_c11, df_c12, df_c21, df_c22],
        }
    else:
        df_bub, df_dew, df_31, df_32 = read_nodsp()
        return {
            "prima_case": prim,
            "has_dsp": dsp,
            "2phase": [df_dew, df_bub],
            "3phase": [df_31, df_32],
        }


def plot_nodsp(stable=False):
    df_bub, df_dew, df_31, df_32 = read_nodsp()
    envels = read_envels()

    prim = envels["prima_case"]
    has_dsp = envels["has_dsp"]

    if stable:
        envels = get_stable(envels, prim, has_dsp)

    df_bub, df_dew = envels["2phase"]
    df_31, df_32 = envels["3phase"]

    index = "T"
    prop = "P"

    show(df_dew, index, prop, label="dew", color="blue")
    show(df_bub, index, prop, label="bub", color="black")

    show(df_31, index, prop, ls="--")
    show(df_32, index, prop, ls="--")


def plot_dsp(stable=True, **kwargs):
    index = "T"
    prop = "P"

    envels = read_envels()

    prim = envels["prima_case"]
    has_dsp = envels["has_dsp"]

    if stable:
        envels = get_stable(envels, has_dsp, prim)

    df_dew, df_bub, df_hpll = envels["2phase"]
    df_c11, df_c12, df_c21, df_c22 = envels["3phase"]

    show(df_dew, index, prop, label="dew", color="blue", **kwargs)
    show(df_bub, index, prop, label="bub", color="black", **kwargs)

    show(df_c11, index, prop, ls="--", **kwargs)
    show(df_c12, index, prop, ls="--", **kwargs)

    show(df_c21, index, prop, ls="--", **kwargs)
    show(df_c22, index, prop, ls="--", **kwargs)

    if prim:
        show(df_hpll, index, prop, label="hpll", color="blue", **kwargs)


def plot(stable=True, **kwargs):
    prim, dsp, three_phase = get_case()
    if LOG: print("Case: ", prim, dsp, three_phase)

    if three_phase:
        if dsp:
            plot_dsp(stable, **kwargs)
        else:
            plot_nodsp(stable)
    else:
        index = "T"
        prop = "P"
        df_dew, df_bub, df_hpll = read_2phase(prim)
        show(df_dew, index, prop, color="blue", **kwargs)
        show(df_bub, index, prop, color="black", **kwargs)
        if prim:
            show(df_hpll, index, prop, color="red", **kwargs)