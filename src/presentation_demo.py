# %%
import matplotlib.ticker as plticker
from matplotlib import pyplot as plt
from matplotlib import rcParams, cycler
from files.map_reading import readReference
from processing.correlate import crossCorrelate

rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c", "k", "c"])


def reference_query_correlation_plot(reference, query, filename, fig1=plt.figure(constrained_layout=True, figsize=(25/2.54, 13/2.54)), suptitle=None):
    gridspec = fig1.add_gridspec(3, 3)

    loc = plticker.MultipleLocator(base=1.0)
    ax1_0 = fig1.add_subplot(gridspec[0, :])
    ax1_1 = fig1.add_subplot(gridspec[1, 0])
    ax1_2 = fig1.add_subplot(gridspec[2, :])

    ax1_0.set_title('Reference')
    ax1_0.bar(range(len(reference)), reference)
    ax1_0.xaxis.set_major_locator(loc)
    ax1_0.set_yticklabels([])
    ax1_0.set_xlim(0, len(reference) - 1)

    ax1_1.set_title('Query')
    ax1_1.bar(range(len(query)), query)
    ax1_1.xaxis.set_major_locator(loc)
    ax1_1.set_yticklabels([])
    ax1_1.set_xlim(0, len(query)-1)

    c = crossCorrelate(reference, query)
    ax1_2.set_title('Cross-correlation')
    ax1_2.xaxis.set_major_locator(loc)
    ax1_2.plot(range(len(c)), c)
    ax1_2.set_xlim(0, len(c)-1)

    if suptitle:
        fig1.suptitle(suptitle)

    fig1.savefig(f"../plots/presentation/{filename}.png",
                 bbox_inches='tight', pad_inches=0)


# %% example
a = [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0]
b = [0, 1, 0, 0, 1, 0, 0, 1, 0]
reference_query_correlation_plot(a, b, "1")

# %%
referenceFile = "../data/ecoli_ref.cmap"
referenceMap = readReference(referenceFile, 1, 0)
fig2, axes2 = plt.subplots(2, 1, constrained_layout=True, figsize=(25/2.54, 13/2.54))


# %% zoomed
fig.axes[0].set_xlim(2400000, 3700000)
fig.set_size_inches((25/2.54, 6/2.54))
fig.savefig("../plots/presentation/plot_molecule1_res10_blur30_errors_True_normalize_False_zoom.png",
            bbox_inches='tight', pad_inches=0)

# %% bluring


def bluring_fig():
    return plt.figure(constrained_layout=True, figsize=(50/2.54, 10/2.54))


before_bluring_data_reference = [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
before_bluring_data_query = [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]

reference_query_correlation_plot(before_bluring_data_reference, before_bluring_data_query, "before_bluring", bluring_fig(), "Before bluring")

after_bluring_data_reference = [1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
after_bluring_data_query = [1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1]

reference_query_correlation_plot(after_bluring_data_reference, after_bluring_data_query, "after_bluring", bluring_fig(), "After bluring")

# %%
