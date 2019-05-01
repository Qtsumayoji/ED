import numpy as np
import pylab as plt
import pandas as pd
import os
import sys
import re
from matplotlib import rc
import matplotlib.cm as cm
import matplotlib.ticker as tick
import seaborn as sns

rc("text", usetex=True)

def RIXS_plot(filename):
    df = pd.read_csv(filename)
    df_pivot = pd.pivot_table(data=df, values="imG", columns="omega", index="k")

    dx = df_pivot.index[1] - df_pivot.index[0]
    dy = df_pivot.columns[1] - df_pivot.columns[0]

    dim_x = len(df_pivot.index)
    dim_y = len(df_pivot.columns)

    x = np.zeros(dim_x + 1)
    y = np.zeros(dim_y + 1)

    x[0] = df_pivot.index[0] - (df_pivot.index[1] - df_pivot.index[0])/2.0
    for i in range(0, dim_x - 1):
        x[i + 1] = df_pivot.index[i] + (df_pivot.index[i + 1] - df_pivot.index[i])/2.0
    x[-1] = df_pivot.index[-1] + (df_pivot.index[-1] - df_pivot.index[-2])/2.0

    y[0] = df_pivot.columns[0] - (df_pivot.columns[1] - df_pivot.columns[0])/2.0
    for i in range(0, dim_y - 1):
        y[i + 1] = df_pivot.columns[i] + (df_pivot.columns[i + 1] - df_pivot.columns[i])/2.0
    y[-1] = df_pivot.columns[-1] + (df_pivot.columns[-1] - df_pivot.columns[-2])/2.0

    fig, ax = plt.subplots(figsize=(8, 6), tight_layout=True)

    fsize = 24
    cmap = "viridis"
    im = ax.pcolormesh(x, y, df_pivot.T, cmap=cmap)
    cb = fig.colorbar(im)
    cb.set_label("Intensity", size=fsize)
    ax.set_xlabel("$k$", size=fsize)
    ax.set_ylabel("$\omega$", size=fsize)

    split = re.split("/", filename)
    dir_result = split[0]

    split = re.split(".csv", split[1])
    fn = split[0]
    plt.savefig(dir_result + "/" + fn + ".eps")

if __name__ == '__main__':
    args = sys.argv
    nfile = len(args) - 1
    
    for i in range(1, nfile + 1):
        filename = args[i]
        RIXS_plot(filename)

    plt.show()
