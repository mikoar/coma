# %%
import matplotlib.ticker as plticker
from matplotlib import cycler
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c", "k", "c"])


def resolution_plot(res1, res10, fig1=plt.figure(constrained_layout=True, figsize=(25 / 2.54, 13 / 2.54)),
                    suptitle=None):

    gridspec = fig1.add_gridspec(3, 2)
    ax1_0 = fig1.add_subplot(gridspec[0, :])
    ax1_1 = fig1.add_subplot(gridspec[1, :])

    loc1 = plticker.MultipleLocator(base=1.0)
    ax1_0.set_title('Resolution = 1')
    ax1_0.bar(range(len(res1)), res1, width=0.95)
    ax1_0.xaxis.set_major_locator(loc1)
    ax1_0.set_yticklabels([])
    ax1_0.set_xlim(-0.5, len(res1) - 0.5)

    loc2 = plticker.MultipleLocator(base=1.0)
    ax1_1.set_title('Resolution = 10')
    ax1_1.bar(range(len(res10)), res10, width=0.95)
    ax1_1.xaxis.set_major_locator(loc2)
    ax1_1.set_yticklabels([])
    ax1_1.set_xlim(-0.5, len(res10) - 0.5)

    if suptitle:
        fig1.suptitle(suptitle)

    fig1.savefig(f"./output/resolution.png",
                 bbox_inches='tight', pad_inches=0)


def res_fig():
    return plt.figure(constrained_layout=True, figsize=(50 / 2.54, 10 / 2.54))


res1 = [1, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
        0, 0, 0, 1, 0, 1, 0, 0, 0, 0]
res10 = [1, 0, 0, 1, 1]

resolution_plot(res1, res10, res_fig())
