function [ExternalStress] = calcExternalField(x, y, z, a, b, c, Nu, EigenStrain, Stiffness)
%   calcExternalField.m 
%   
%   Calculate elastic stress outside the ellipsoidal inclusion 
%
%   Arguments:
%       x, y, z = coordinates 
%       a, b, c = semi-axes of ellipsoid
%       Nu = Poisson's ratio of matrix 
%       EigenStrain = total eigen strain of inclusion, in Voigt order 
%
%   Based on the formulations originally derived by Ju & Sun:
%       Ju & Sun, 1999. Journal of Applied Mechanics, 66, p570. 
%       Ju & Sun, 2001. Intl Journal of Solids & Structures, 38, p183.  
%
%   David Healy
%   April 2008 
%
%   Please let me know about any bugs or errors: 
%       d.healy@curtin.edu.au
%       david.healy@mac.com 
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.

%   calculate Eshelby G tensor for this point 
Eshelby_G = calcEshelbyG(x, y, z, a, b, c, Nu) ; 

ExternalStress = ( Stiffness .* Eshelby_G ) * EigenStrain ; 