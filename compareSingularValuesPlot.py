#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import inspect, os
from scipy.io import netcdf
import sys

myname = inspect.getfile(inspect.currentframe())
print("This is " + myname)
print("Usage:")
print("  " + myname + " <List of 1 or more bdistrib_out.XXX.nc files>")
print()

if len(sys.argv) < 2:
    print(
        "Error! You must list at least one bdistrib_out.XXX.nc file as a command-line argument."
    )
    exit(1)


def maximizeWindow():
    # Maximize window. The command for this depends on the backend.
    mng = plt.get_current_fig_manager()
    try:
        mng.resize(*mng.window.maxsize())
    except AttributeError:
        try:
            mng.window.showMaximized()
        except AttributeError:
            pass


svd_s_transferMatrix_many = []
svd_s_inductance_plasma_middle_many = []
dataNames = []

for whichFile in range(1, len(sys.argv)):
    filename = sys.argv[whichFile]
    f = netcdf.netcdf_file(filename, "r", mmap=False)
    svd_s_transferMatrix = f.variables["svd_s_transferMatrix"][()]
    svd_s_inductance_plasma_middle = f.variables["svd_s_inductance_plasma_middle"][()]
    pseudoinverse_thresholds = f.variables["pseudoinverse_thresholds"][()]
    f.close()

    print("Read data from file " + filename)

    svd_s_inductance_plasma_middle_many.append(svd_s_inductance_plasma_middle)

    for whichThreshold in range(len(pseudoinverse_thresholds)):
        svd_s_transferMatrix_many.append(svd_s_transferMatrix[whichThreshold, :])
        dataNames.append(
            filename + " (thresh=" + str(pseudoinverse_thresholds[whichThreshold]) + ")"
        )


##########################################################
# Make plot for inductance matrix
##########################################################

fig = plt.figure(1)
fig.patch.set_facecolor("white")

for whichFile in range(1, len(sys.argv)):
    filename = sys.argv[whichFile]
    plt.semilogy(
        svd_s_inductance_plasma_middle_many[whichFile - 1], ".-", label=filename
    )
plt.title("Singular values of plasma-control inductance matrix")
plt.legend(frameon=False, prop=dict(size="x-small"), loc=3)
plt.grid(True)

titleString = (
    "Plot generated by "
    + os.path.abspath(inspect.getfile(inspect.currentframe()))
    + "\nRun in "
    + os.getcwd()
)
ax = fig.add_axes((0, 0, 1, 1), frameon=False)
ax.text(0.5, 0.99, titleString, horizontalalignment="center", verticalalignment="top")

maximizeWindow()

##########################################################
# Make plot for transfer matrix
##########################################################

fig = plt.figure(2)
fig.patch.set_facecolor("white")

for i in range(len(dataNames)):
    plt.semilogy(svd_s_transferMatrix_many[i], ".-", label=dataNames[i])
plt.title("Singular values of transfer matrix")
plt.legend(frameon=False, prop=dict(size="x-small"), loc=3)
plt.grid(True)

titleString = (
    "Plot generated by "
    + os.path.abspath(inspect.getfile(inspect.currentframe()))
    + "\nRun in "
    + os.getcwd()
)
ax = fig.add_axes((0, 0, 1, 1), frameon=False)
ax.text(0.5, 0.99, titleString, horizontalalignment="center", verticalalignment="top")

maximizeWindow()

plt.show()
