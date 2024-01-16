The SSG Reynolds stress model is apprently suitable to obtain 
the log-layer solution without wall damping functions of the pressure-strain 
tensor (Wilcox 2006). However, it works only near the no slip solid 
boundaries with steep velocity gradients. 
The available SSG model in OpenFOAM does not include the non-linear free surface 
damping functions and dissipation rate boundary condition, which are needed
to simulate turbulence-driven secondary currents. 
Therefore, these functions were added in the modified SSG model or SSGMod based 
on the literature and validated for a supercritical uniform case with 
rectangular cross-section before studying the effect of numerous cross-sectional 
geometries on narrow channel flows comparable to sediment bypass tunnels.  

The modifications were tested in -Dev version. They should be working 
in OpenFOAM-10 also. 

For more details follow the the following manuscript: 
Kadia, S., Lia, L., Albayrak, I., and Pummer, E. (2024). “The effect of cross-sectional geometry 
on the high-speed narrow open channel flows: An updated Reynolds stress model study.” Computers & Fluids, in press.

References:
1. Wilcox, D. C. (2006). Turbulence Modelling for CFD. DCW Industries, Inc., La Canada, California.
