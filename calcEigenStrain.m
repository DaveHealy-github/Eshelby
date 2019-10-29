function [EigenStrain] = calcEigenStrain(M_C, I_C, S, M_Strain, I_Strain)
%   calcEigenStrain.m
%
%   Calculate inclusion eigenstrain 
%
%   Arguments:
%       M_C = matrix stiffness, 6 x 6 matrix
%       I_C = inclusion stiffness, 6 x 6 matrix
%       S = Eshelby S tensor, 6 x 6 matrix 
%       M_Strain = initial matrix strain, 6 x 1 vector in Voigt order  
%       I_Strain = initial inclusion strain, 6 x 1 vector in Voigt order  
%
%   Based on formulation in: 
%       Ju & Sun, 1999. Journal of Applied Mechanics, 66, p570. 
%       Ju & Sun, 2001. Intl Journal of Solids & Structures, 38, p183.  
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

%   calculate difference between inclusion and matrix stiffness matrices 
%   Ju & Sun 1999, eqns 7-9
Delta_C = I_C - M_C ; 

%   get the inverse of this difference 
Inv_Delta_C = inv(Delta_C) ; 

%   multiply by matrix stiffness for matrix A
A = Inv_Delta_C .* M_C ; 

%   and by inclusion stiffness for matrix B
B = Inv_Delta_C .* I_C ; 

%   calculate S + A
Sum = S + A ; 

%   find the inverse 
Inv_Sum = inv(Sum) ; 

%   check if prescribed inclusion strain is 0 
if ( I_Strain == 0 ) 
    disp('I_Strain is 0') ; 
    %   Ju & Sun 1999, eqn 9b (in the text)
    EigenStrain = (-Inv_Sum) * M_Strain ; 
else 
    disp('I_Strain is not 0') ; 
    %   Ju & Sun 1999, eqn 7  
    EigenStrain = Inv_Sum * ( ( B * I_Strain ) - M_Strain ) ; 
end ; 