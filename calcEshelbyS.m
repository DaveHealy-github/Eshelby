function [S] = calcEshelbyS(Nu, a, b, c)
%   calcEshelbyS.m 
%
%   Calculates Eshelby S tensor components in a reduced 6x6 matrix. 
%
%   Based on the formulation in:
%       Ju & Sun, 2001. International Journal of Solids and Structures, 38, p183. 
%
%   Arguments:
%       Nu = Poisson's ratio of matrix
%       a, b, c = semi-axes of ellipsoid 
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

disp(' ') ; 
disp('Calculating Eshelby S tensor...') ; 

S = zeros(6,6) ; 
ArrayS = zeros(3,3,3,3) ;  

%   some useful expressions... 
KroneckerDelta = [ 1 0 0 ; ... 
                   0 1 0 ; ... 
                   0 0 1 ] ; 
Nu2 = Nu * 2 ; 
Nu4 = Nu * 4 ; 
Nu5 = Nu * 5 ; 
OneMinusNu = 1 - Nu ; 
Alpha = a / c ; 
AlphaSquared = Alpha^2 ; 
AlphaSquaredMinusOne = AlphaSquared - 1 ; 
OneMinusAlphaSquared = 1 - AlphaSquared ; 

%   S defined differently for different ellipsoid shapes...
if ( Alpha ~= 1 )  

    if ( Alpha > 1 )   
        
        %   prolate spheroid: Alpha > 1 
        disp('Prolate spheroidal inclusion...') ; 
        
        %   Ju & Sun 2001, eqn A.19a 
        g0 = ( Alpha / AlphaSquaredMinusOne^1.5 ) ...
                * ( log( Alpha + sqrt(AlphaSquaredMinusOne) ) ...
                    - ( Alpha * sqrt(AlphaSquaredMinusOne) ) ) ; 
    
    else

        %   oblate spheroid: Alpha < 1
        disp('Oblate spheroidal inclusion...') ;
        
        %   Ju & Sun 2001, eqn A.19b 
        g0 = ( Alpha / OneMinusAlphaSquared^1.5 ) ...
                * ( ( Alpha * sqrt(OneMinusAlphaSquared) ) ...
                   - ( acos(Alpha) ) ) ;   

    end ;
    
    %   S1 components - Ju & Sun 2001, eqns A.12-A.15 
    S1_11 = ( ( Nu4 + ( 2 / AlphaSquaredMinusOne ) ) * g0 ) ...
                + Nu4 + ( 4 / ( 3 * AlphaSquaredMinusOne ) ) ; 
    S1_12 = ( ( Nu4 - ( ( 2 * AlphaSquared + 1 ) / AlphaSquaredMinusOne ) ) * g0 ) ...
                + Nu4 - ( ( 2 * AlphaSquared ) / AlphaSquaredMinusOne ) ; 
    S1_13 = S1_12 ;      

    S1_21 = ( ( -Nu2 - ( ( 2 * AlphaSquared + 1 ) / AlphaSquaredMinusOne ) ) * g0 ) ...
                - ( ( 2 * AlphaSquared ) / AlphaSquaredMinusOne ) ; 
    S1_22 = ( ( -Nu2 + ( ( 4 * AlphaSquared - 1 ) / ( 4 * AlphaSquaredMinusOne ) ) ) * g0 ) ...
                + ( AlphaSquared / ( 2 * AlphaSquaredMinusOne ) ) ; 
    S1_23 = S1_22 ; 
    
    S1_31 = S1_21 ; 
    S1_32 = S1_22 ; 
    S1_33 = S1_22 ; 
    
    %   S2 components - Ju & Sun 2001, eqns A.16-A.18
    S2_11 = ( ( -Nu4 + ( ( 4 * AlphaSquared - 2 ) / AlphaSquaredMinusOne ) ) * g0 ) ...
                - Nu4 + ( ( 12 * AlphaSquared - 8 ) / ( 3 * AlphaSquaredMinusOne ) ) ; 
    S2_12 = ( ( -Nu - ( ( AlphaSquared + 2 ) / AlphaSquaredMinusOne ) ) * g0 ) ...
                - Nu2 - ( 2 / AlphaSquaredMinusOne ) ; 
    S2_13 = S2_12 ; 
    
    S2_21 = S2_12 ; 
    S2_22 = ( ( Nu2 - ( ( 4 * AlphaSquared - 7 ) / ( 4 * AlphaSquaredMinusOne ) ) ) * g0 ) ...
                + ( AlphaSquared / ( 2 * AlphaSquaredMinusOne ) ) ; 
    S2_23 = S2_22 ; 
    
    S2_31 = S2_12 ; 
    S2_32 = S2_22 ; 
    S2_33 = S2_22 ; 
    
    %   map these elements into two 3x3 arrays 
    S1 = [ S1_11 S1_12 S1_13 ; ...
           S1_21 S1_22 S1_23 ; ...
           S1_31 S1_32 S1_33 ] ; 
            
    S2 = [ S2_11 S2_12 S2_13 ; ...
           S2_21 S2_22 S2_23 ; ...
           S2_31 S2_32 S2_33 ] ; 
    
    %   Ju & Sun 2001, eqn 23
    for i = 1:3 
        for j = 1:3 
            for k = 1:3
                for l = 1:3
                    ArrayS(i,j,k,l) = ( 1 / ( 4 * OneMinusNu ) ) ...
                                      * ( ( S1(i,k) * KroneckerDelta(i,j) * KroneckerDelta(k,l) ) ...
                                        + ( S2(i,j) * ( KroneckerDelta(i,k) * KroneckerDelta(j,l) ...
                                                      + KroneckerDelta(i,l) * KroneckerDelta(j,k) ) ) ) ;
                end ; 
            end ; 
        end ; 
    end ;      
    
else 
    %   sphere: Alpha = 1
    disp('Spherical inclusion...') ; 
    
    %   Ju & Sun 2001, eqn A.20 
    for i = 1:3 
        for j = 1:3 
            for k = 1:3
                for l = 1:3
                    ArrayS(i,j,k,l) = ( 1 / ( 15 * OneMinusNu ) ) ...
                                        * ( ( Nu5 - 1 ) ...
                                          * KroneckerDelta(i,j) * KroneckerDelta(k,l) ...
                                        + ( 4.0 - Nu5 ) ...
                                          * ( KroneckerDelta(i,k) * KroneckerDelta(j,l) ...
                                            + KroneckerDelta(i,l) * KroneckerDelta(j,k) ) ) ; 
    
                end ; 
            end ; 
        end ; 
    end ;      

end ; 

% now map the 4D 3x3x3x3 array into the 2D 6x6 matrix
S(1,1) = ArrayS(1,1,1,1) ; 
S(2,2) = ArrayS(2,2,2,2) ; 
S(3,3) = ArrayS(3,3,3,3) ; 
S(4,4) = ArrayS(2,3,2,3) ; 
S(5,5) = ArrayS(3,1,3,1) ; 
S(6,6) = ArrayS(1,2,1,2) ; 
S(1,2) = ArrayS(1,1,2,2) ; 
S(1,3) = ArrayS(1,1,3,3) ; 
S(2,1) = ArrayS(2,2,1,1) ; 
S(3,1) = ArrayS(3,3,1,1) ; 
S(3,2) = ArrayS(3,3,2,2) ; 
S(2,3) = ArrayS(2,2,3,3) ; 
