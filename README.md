# Surface-Integral-Equation-program
Calculating optical cross sections from an arbitrary scatterer using surface integral equation.

Run the code through runner.m.

How to build geometrical mesh by COMSOL:

Create the 3D object in COMSOL
Create the free triangular mesh at the object`s surface
Export the surface mesh in .mphtxt format
Import the created mesh in MATLAB. Import these variables: p (triangle vertex points), the first chunk of columns consists of three rows, and t (triangle connection between the vertex points), the second of the three row columns. Add variable tby 1. Save these two variables and create a line to load these variables in runner.m
