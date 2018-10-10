# Hybrid Space Plasma plotting

## Simulation
This code will only read and plot the data provided by a simulation.  The simulation this code is based off of can be found here: https://github.com/pdelamere/Hybrid_3d_buf

## Organization
Paramaters are read in using Read_para.pro, must be run before Create_movie

Coordinates are read in and the grid is constructed using Read_coord.pro

Read_scalar and Read_vector read in particle data

Create_movie is the main program, creates each frame and animates it

img_cont is based on the built in IDL procedure image_cont, graphs the frame and outputs it to the display

## Known issues
No color bars are created, and graphs are not labelled.

## Contact
This version of the code written by James.Mahon@lasp.colorado.edu
