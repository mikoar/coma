# %%
import matplotlib.ticker as plticker
from matplotlib import pyplot as plt
from matplotlib import rcParams, cycler
from files.map_reading import readReference
from processing.correlate import crossCorrelate


# %% example
a = [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0]
b = [0, 1, 0, 0, 1, 0, 0, 1, 0]
c = crossCorrelate(a, b)


fig1 = plt.figure(constrained_layout=True, figsize=(25/2.54, 13/2.54))
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c", "k", "c"])
gridspec = fig1.add_gridspec(3, 3)

loc = plticker.MultipleLocator(base=1.0)
ax1_0 = fig1.add_subplot(gridspec[0, :])
ax1_1 = fig1.add_subplot(gridspec[1, 0])
ax1_2 = fig1.add_subplot(gridspec[2, :])

ax1_0.set_title('Reference')
ax1_0.bar(range(len(a)), a)
ax1_0.xaxis.set_major_locator(loc)
ax1_0.set_yticklabels([])
ax1_0.set_xlim(0, len(a) - 1)

ax1_1.set_title('Query')
ax1_1.bar(range(len(b)), b)
ax1_1.xaxis.set_major_locator(loc)
ax1_1.set_yticklabels([])
ax1_1.set_xlim(0, len(b)-1)

ax1_2.set_title('Cross-correlation')
ax1_2.xaxis.set_major_locator(loc)
ax1_2.plot(range(len(c)), c)
ax1_2.set_xlim(0, len(c)-1)

fig1.savefig("../plots/presentation/1.png",
             bbox_inches='tight', pad_inches=0)

# %%
referenceFile = "../data/ecoli_ref.cmap"
referenceMap = readReference(referenceFile, 1, 0)
fig2, axes2 = plt.subplots(2, 1, constrained_layout=True, figsize=(25/2.54, 13/2.54))


# %% zoomed
fig.axes[0].set_xlim(2400000, 3700000)
fig.set_size_inches((25/2.54, 6/2.54))
fig.savefig("../plots/presentation/plot_molecule1_res10_blur30_errors_True_normalize_False_zoom.png",
            bbox_inches='tight', pad_inches=0)
