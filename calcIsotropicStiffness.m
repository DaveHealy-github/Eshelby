function [Stiffness] = calcIsotropicStiffness(E, Nu)
%   calcIsotropicStiffness.m 
%
%   Populate 6 x 6 isotropic elastic stiffness matrix
%
%   Arguments:
%       E = Young's modulus
%       Nu = Poisson's ratio
%
%   David Healy
%   October 2008 
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

Lambda = ( E * Nu ) / ( ( 1 + Nu ) * ( 1 - 2 * Nu ) ) ; 
OnDiagonal = Lambda + ( 2 * E / ( 2 * ( 1 + Nu ) ) ) ; 
OffDiagonal = Lambda ; 

Stiffness = [ OnDiagonal OffDiagonal OffDiagonal zeros(1, 3) ; ...
              OffDiagonal OnDiagonal OffDiagonal zeros(1, 3) ; ...
              OffDiagonal OffDiagonal OnDiagonal zeros(1, 3) ; ...
              zeros(1, 3) ( (OnDiagonal - OffDiagonal) / 2 ) zeros(1, 2) ; ... 
              zeros(1, 4) ( (OnDiagonal - OffDiagonal) / 2 ) 0 ; ...
              zeros(1, 5) ( (OnDiagonal - OffDiagonal) / 2 ) ] ; 
