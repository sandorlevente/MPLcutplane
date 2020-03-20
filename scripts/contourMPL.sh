#!/bin/bash
# vtk to matplotlib plotting script
# run from case/scripts folder!
echo "Plotting script started"
#creating directories
mkdir ../postProcessing/figures
mkdir ../postProcessing/triangulatedCuttingPlanes
echo "Directories created"

#activating paraview python environment
pvpython triangulateCutplanes.py
echo "Triangulation finished"

#activating the anaconda environment with matplotlib (only on linux)
source ~/anaconda3/etc/profile.d/conda.sh
conda activate base

python plotPlanes.py

echo "Plotting script finished"