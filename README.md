# Orbital Dynamics Processes
-------------------

This is a python module that contains a few useful calculations an conversions of 
useful parameters in orbital dynamics. It contains the following (complete) functions:

* `el2xv` -- Function to change orbital elements semimajor axis (AU), eccentricity, inclination (rad), longitude of ascending node - capom (rad), argument of pericenter - omega (rad), and the mean anomaly - capm (rad) to position (AU) and velocity (AU/day) using mu (AU^3/day^2 system). Only tested for the elliptical case.
* `xv2el` -- Function to change the state vectors (position and velocity) to the orbital elements [semimajor axis, eccentricity,  inclination, longitude of ascending node (capom), argument of pericenter (omega), and the mean anomaly (capm)] Only tested for the elliptical case
* `orb_dist_Dd` -- Function to calculate the Drummond (1981) Orbital Similarity distance (Dd) given 6 orbital elements of two bodies:
    Semimajor axis (AU), Eccentricity, Inclination (rad), Longitude of Ascending Node (rad), & Argument of periapsis (rad), & Mean Anomaly (rad) (transferred in tuples) and the mass of the sun. Equations were taken from Rozek et al. (2010)
