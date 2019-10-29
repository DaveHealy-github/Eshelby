*********************************************************
readme.txt for Eshelby's solution in MATLAB

David Healy
Curtin University
Perth, Australia

d.healy@curtin.edu.au 
david.healy@mac.com

October 2008 
*********************************************************

Please note the following conventions and points about the code:

1. Ellipsoid semi-axes a, b and c are parallel to coordinate axes x, y and z respectively. 

2. Ellipsoid is centred at (0, 0, 0). 

3. The published exact solutions for Eshelby's formulation *only* apply to spheroids i.e. spheres, oblate spheroids or prolate spheroids. 

4. For non-spherical ellipsoids, semi-axis a *must* be either the shortest or the longest semi-axis. The ratio a/c is used in the classification of the ellipsoid when computing the S and G tensors.

5. For non-spheroidal ellipsoids, the solution for the elastic field is not in a closed form, i.e. not expressed through analytical expressions, but depends on inexact elliptic integrals. At the time of writing, these elliptic integrals are not a part of the MATLAB function library.

6. Tensile stress and stretching strain are positive, compressive stress and shortening strain are negative.

7. If you find any bugs/errors, please let me know: d.healy@curtin.edu.au

