import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl

fig, ax = pl.subplots(2, 2, sharex=True, sharey="row")

file = h5py.File("test000.hdf5", "r")
box = file["/Header"].attrs["BoxSize"]
coords = file["/PartType0/Coordinates"]
dens = file["/PartType0/TestDensity"]
dsum = file["/PartType0/TestNeighbourDensitySum"]

r = np.sqrt(
    (coords[:, 0] - box[0]) ** 2
    + (coords[:, 1] - box[1]) ** 2
    + (coords[:, 2] - box[2]) ** 2
)

ax[0][0].plot(r, dens, ".")
ax[1][0].plot(r, dsum, ".")

ax[0][0].set_title("Serial loop")

file = h5py.File("test001.hdf5", "r")
box = file["/Header"].attrs["BoxSize"]
coords = file["/PartType0/Coordinates"]
dens = file["/PartType0/TestDensity"]
dsum = file["/PartType0/TestNeighbourDensitySum"]

r = np.sqrt(
    (coords[:, 0] - box[0]) ** 2
    + (coords[:, 1] - box[1]) ** 2
    + (coords[:, 2] - box[2]) ** 2
)

ax[0][1].plot(r, dens, ".")
ax[1][1].plot(r, dsum, ".")

ax[0][1].set_title("Task-based")

ax[0][0].set_ylabel("density copy")
ax[1][0].set_ylabel("neighbour sum")

ax[1][0].set_xlabel("r (m)")
ax[1][1].set_xlabel("r (m)")

pl.tight_layout()
pl.savefig("results.png", dpi=300)
