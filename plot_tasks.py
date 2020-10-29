################################################################################
# This file is part of CMacIonize
# Copyright (C) 2018, 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# CMacIonize is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CMacIonize is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with CMacIonize. If not, see <http://www.gnu.org/licenses/>.
################################################################################

##
# @file plot_tasks.py
#
# @brief Script to plot the task plot for a given file with task output.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# import modules
import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib.ticker import AutoMinorLocator
import pylab as pl
import sys
import argparse

# parse the command line arguments
argparser = argparse.ArgumentParser(
    description="Plot task plot based on a given task output file."
)

argparser.add_argument("-n", "--name", action="store", required=True)

args = argparser.parse_args(sys.argv[1:])

name = args.name

# change the default matplotlib settings to get nicer plots
pl.rcParams["figure.figsize"] = (12, 10)
pl.rcParams["font.size"] = 14

# name labels and colours for the various task types
# add extra colours and labels here if new tasks are created
task_names = [
    "test density",
    "test neighbour density sum internal",
    "test neighbour density sum external",
    "test neighbour density sum ghost",
]
task_colors = pl.cm.ScalarMappable(cmap="tab20").to_rgba(
    np.linspace(0.0, 1.0, len(task_names))
)

# load the data
print("Plotting tasks for", name, "...")
data = np.loadtxt(
    name,
    dtype={
        "names": ("index", "type", "thread", "start", "stop"),
        "formats": ("u4", "u1", "i4", "u8", "u8"),
    },
)
# tasks 0-18 are not used, so we rebase the numbers here
data["type"] -= 18

task_flags = [
    len(data[data["type"] == task]) > 0 for task in range(len(task_names))
]

# get information about the system
nthread = int(data["thread"].max()) + 1

# get the minimum and maximum time stamp and compute the time to fraction
# conversion factor for each node
tmin = data["start"].min()
tmax = data["stop"].max()

## make the plot

fig, ax = pl.subplots(1, 1, sharex=True)

ax.axvline(x=0, linestyle="--", color="k", linewidth=0.8)
ax.axvline(x=tmax - tmin, linestyle="--", color="k", linewidth=0.8)

# now plot the tasks
alltime = 0
# loop over the threads
for i in range(nthread):
    # filter out the data for this thread
    thread = data[data["thread"] == i]

    # create the task plot
    bar = [
        ((task["start"] - tmin), (task["stop"] - task["start"]))
        for task in thread
    ]
    colors = [task_colors[task["type"]] for task in thread]
    ax.broken_barh(bar, (i - 0.4, 0.8), facecolors=colors, edgecolor="none")

# add empty blocks for the legend
for i in range(len(task_colors)):
    if task_flags[i]:
        ax.plot([], [], color=task_colors[i], label=task_names[i])

# add the legend and clean up the axes
ax.legend(loc="upper center", ncol=min(2, len(task_colors)))
ax.set_ylim(-1.0, nthread * 1.1)
ax.set_yticks([])

ax.set_xlabel("Simulation time (CPU ticks)")

ax.xaxis.set_minor_locator(AutoMinorLocator())

# finalize and save the plot
pl.tight_layout()
pl.savefig("{0}.png".format(name[:-4]))
