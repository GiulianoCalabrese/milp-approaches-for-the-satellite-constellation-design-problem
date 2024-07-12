# MILP Approaches for the Satellite Constellation Design Problem

Code giving results from the paper with the same title :

Mencarelli, L., Floquet, J., Georges, F. et al. Mixed integer (non)linear approaches for the satellite constellation design problem. Optim Eng (2022). https://doi.org/10.1007/s11081-022-09774-9

You need install library JuMP, Cbc and SatelliteToolbox
You need launch with Julia the Keplerian file and for validate the results you use the Validation code.

## Getting started

Abstract In this paper, we propose mathematical optimization models to solve the satellite constellation design problem for discontinuous coverage. In such a design problem, the aim is to determine the minimal number of satellites (and, incidentally, their 3D placements) in order to observe a fixed Earth region within a given revisiting time. Two Mixed Integer Nonlinear formulations are introduced. The first one is a feasibility problem based on the direct mathematical definition of pixel observability. The second one consists in introducing a set of indicator variables which specify if a satellite observes a
pixel at a given time-stamp. In order to obtain a linear problem, the possible positions of the satellites are discretized. Finally, computational results show the potential and limitations of the proposed approaches.

## Authors and acknowledgment
Luca Mencarelli, Julien Floquet, Frederic Georges

