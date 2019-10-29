function [EshelbyG] = calcEshelbyG(x, y, z, a, b, c, Nu)
%   calcEshelbyG.m
%   
%   Calculates Eshelby G tensor components in a reduced 6 x 6 matrix
%
%   Based on the formulation in:
%       Ju & Sun, 2001. International Journal of Solids and Structures, 38, p183. 
%
%   *** with corrections to equations A.1 and A.2 
%
%   Arguments:
%       x, y, z = coordinate position of a point outside the inclusion
%       a, b, c = semi-axes of ellipsoidal inclusion 
%       Nu = Poisson's ratio of the matrix 
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

ArrayGbar = zeros(3, 3, 3, 3) ; 
ArrayS7 = zeros(3, 3, 3, 3) ; 

%   some useful, but tedious, expressions... 
KroneckerDelta = [ 1 0 0 ; ... 
                   0 1 0 ; ... 
                   0 0 1 ] ; 
Nu2 = Nu * 2 ; 
Nu4 = Nu * 4 ; 
OneMinusNu = 1 - Nu ; 

Alpha = a / c ; 
AlphaSquared = Alpha^2 ; 
AlphaSquaredMinusOne = AlphaSquared - 1 ; 
OneMinusAlphaSquared = 1 - AlphaSquared ; 

xSquared = x^2 ; 
ySquared = y^2 ; 
zSquared = z^2 ; 
rSquared = xSquared + ySquared + zSquared ; 

aSquared = a^2 ; 
bSquared = b^2 ; 
cSquared = c^2 ; 

Lamda = ( rSquared - aSquared - bSquared ... 
            + sqrt( ( ( ( rSquared + aSquared - bSquared ) ...
            * ( rSquared + aSquared - bSquared ) ) ...
            - ( 4 * xSquared * ( aSquared - bSquared ) ) ) ) ) ...
            / 2 ; 

Theta = [ x / ( aSquared + Lamda ) ; ...
          y / ( bSquared + Lamda ) ; ...
          z / ( cSquared + Lamda ) ] ; 
      
ThetaAll = Theta(1)^2 + Theta(2)^2 + Theta(3)^2 ; 

Rho = [ a / sqrt( a^2 + Lamda ) ; ...
        b / sqrt( b^2 + Lamda ) ; ...
        c / sqrt( c^2 + Lamda ) ] ; 
    
Nhat = [ x / ( ( a^2 + Lamda ) * sqrt(ThetaAll) ) ; ...
         y / ( ( b^2 + Lamda ) * sqrt(ThetaAll) ) ; ...
         z / ( ( c^2 + Lamda ) * sqrt(ThetaAll) ) ] ; 
     
RhoAll = Rho(1)^2 + Rho(2)^2 + Rho(3)^2 ; 
RhoCubed = Rho(1) * Rho(2) * Rho(3) ; 
RhoSquaredThetaAll = ( Rho(1)^2 * Theta(1)^2 ) ... 
                   + ( Rho(2)^2 * Theta(2)^2 ) ... 
                   + ( Rho(3)^2 * Theta(3)^2 ) ;
               
%   Eshelby G tensor depends on ellipsoid shape...
if ( Alpha ~= 1 ) 
    
    if ( Alpha > 1 ) 
        
        %   prolate spheroid 
        %disp('Prolate spheroidal inclusion...') ; 
        
        %   Ju & Sun 2001, eqn A.8a 
        g = ( ( -AlphaSquared / AlphaSquaredMinusOne ) ... 
             * ( Rho(2)^2 ) / Rho(1) ) ...
             + ( ( Alpha / AlphaSquaredMinusOne^1.5 ) ...
             * log( ( sqrt(AlphaSquaredMinusOne) * Rho(2) ) ...
             + ( ( Alpha * Rho(2) ) / Rho(1) ) ) ) ; 
         
    else 
        
        %   oblate spheroid 
        %disp('Oblate spheroidal inclusion...') ; 
        
        %   Ju & Sun 2001, eqn A.8b
        g = ( ( AlphaSquared / OneMinusAlphaSquared ) ...
             * ( Rho(2)^2 ) / Rho(1) ) ...
             - ( ( Alpha / OneMinusAlphaSquared^1.5 ) ...
             * ( acos( Alpha * Rho(2) / Rho(1) ) ) ) ; 
        
    end ; 

    %   S components of Gbar 
    %   Ju & Sun 2001, eqns A.1 - A.7
    %   NB: A.1 corrected - reverse signs in first term on RHS
    %   NB: A.2 corrected - reverse signs in first term on RHS
    %   NB: - these corrections are printed incorrectly (!) in 
    %       Healy et al., 2006, EPSL. 
    S1_11 = ( ( Nu4 + ( 2 / AlphaSquaredMinusOne ) ) * g ) ... 
              - ( ( 2 / ( 3 * AlphaSquaredMinusOne ) ) * Rho(1)^3 ) ... 
              + ( ( Nu4 + 2 / AlphaSquaredMinusOne ) * Rho(1) * Rho(2)^2 ) ;
    S1_12 = ( ( Nu4 - ( ( 2 * AlphaSquared + 1 ) / AlphaSquaredMinusOne ) ) * g ) ...
              + ( ( Nu4 - ( ( 2 * AlphaSquared ) / AlphaSquaredMinusOne ) ) ...
              * Rho(1) * Rho(2)^2 ) ; 
    S1_13 = S1_12 ; 
    
    S1_21 = ( ( -Nu2 - ( ( 2 * AlphaSquared + 1 ) / AlphaSquaredMinusOne ) ) * g ) ...
              - ( ( ( 2 * AlphaSquared ) / AlphaSquaredMinusOne ) ...
              * Rho(1) * Rho(2)^2 ) ; 
    S1_22 = ( ( -Nu2 + ( ( 4 * AlphaSquared - 1 ) ... 
                / ( 4 * AlphaSquaredMinusOne ) ) ) * g ) ... 
                + ( ( AlphaSquared / ( 2 * AlphaSquaredMinusOne ) ) ... 
                * ( Rho(2)^4 ) / Rho(1) ) ; 
    S1_23 = S1_22 ; 
    
    S1_31 = S1_21 ; 
    S1_32 = S1_22 ; 
    S1_33 = S1_22 ; 
    
    S2_11 = ( ( -Nu4 + ( ( 4 * AlphaSquared - 2 ) / AlphaSquaredMinusOne ) ) * g ) ...
                - ( ( 2 / ( 3 * AlphaSquaredMinusOne ) ) * Rho(1)^3 ) ...
                - ( ( Nu4 - ( 4 * AlphaSquared - 2 ) / AlphaSquaredMinusOne ) ...
                * Rho(1) * Rho(2)^2 ) ; 
    S2_12 = ( ( -Nu - ( ( AlphaSquared + 2 ) / AlphaSquaredMinusOne ) ) * g ) ...
                - ( ( Nu2 + 2 / AlphaSquaredMinusOne ) ... 
                * Rho(1) * Rho(2)^2 ) ; 
    S2_13 = S2_12 ; 
    
    S2_21 = S2_12 ; 
    S2_22 = ( ( Nu2 - ( ( 4 * AlphaSquared - 7 ) ...
                / ( 4 * AlphaSquaredMinusOne ) ) ) * g ) ... 
                + ( ( AlphaSquared / ( 2 * AlphaSquaredMinusOne ) ) ... 
                * Rho(2)^4 / Rho(1) ) ; 
    S2_23 = S2_22 ; 
    
    S2_31 = S2_12 ; 
    S2_32 = S2_22 ; 
    S2_33 = S2_22 ; 
    
    %   map S1 and S2 components into 3 x 3 matrices  
    S1 = [ S1_11 S1_12 S1_13 ; ...
           S1_21 S1_22 S1_23 ; ... 
           S1_31 S1_32 S1_33 ] ; 
    S2 = [ S2_11 S2_12 S2_13 ; ... 
           S2_21 S2_22 S2_23 ; ...
           S2_31 S2_32 S2_33 ] ; 
    
    %   S3 components, Ju & Sun 2001, eqn 12 
    S3 = [ ( 2 * RhoCubed ) * ( 1 - Rho(1)^2 ) ; ... 
           ( 2 * RhoCubed ) * ( 1 - Rho(2)^2 ) ; ...
           ( 2 * RhoCubed ) * ( 1 - Rho(3)^2 ) ] ; 

    %   S4 components, Ju & Sun 2001, eqn 13 
    S4 = [ ( 2 * RhoCubed ) * ( 1 - Nu2 - Rho(1)^2 ) ; ... 
           ( 2 * RhoCubed ) * ( 1 - Nu2 - Rho(2)^2 ) ; ...
           ( 2 * RhoCubed ) * ( 1 - Nu2 - Rho(3)^2 ) ] ; 
    
    %   S5 components, Ju & Sun 2001, eqn 14 
    S5 = [ ( 2 * RhoCubed ) * ( Nu - Rho(1)^2 ) ; ... 
           ( 2 * RhoCubed ) * ( Nu - Rho(2)^2 ) ; ...
           ( 2 * RhoCubed ) * ( Nu - Rho(3)^2 ) ] ; 
    
    %   S6 components, Ju & Sun 2001, eqn 15
    S6 = [ ( 2 * RhoCubed ) * ( Nu - Rho(1)^2 ) ; ... 
           ( 2 * RhoCubed ) * ( Nu - Rho(2)^2 ) ; ...
           ( 2 * RhoCubed ) * ( Nu - Rho(3)^2 ) ] ; 
    
    %   81 elements of S7, Ju & Sun 2001, eqn 16 
    for i = 1:3 
        for j = 1:3 
            for k = 1:3 
                for l = 1:3 
                    ArrayS7(i,j,k,l) = ( 2 * RhoCubed ) ...
                            * ( 2 * ( Rho(i)^2 + Rho(j)^2 + Rho(k)^2 + Rho(l)^2 ) ...
                            + RhoAll - ( ( 4 * RhoSquaredThetaAll ) / ThetaAll ) ...
                            - 5 ) ; 
                end ; 
            end ; 
        end ; 
    end ; 
    
    %   81 elements of Gbar, Ju & Sun 2001, eqn 9 
    for i = 1:3 
        for j = 1:3 
            for k = 1:3 
                for l = 1:3 
                    ArrayGbar(i,j,k,l) = ( 1 / ( 4 * OneMinusNu ) ) ... 
                             * ( ( S1(i,k) * KroneckerDelta(i,j) * KroneckerDelta(k,l) ) ... 
                               + ( S2(i,j) * ( ( KroneckerDelta(i,k) * KroneckerDelta(j,l) ) ...
                                             + ( KroneckerDelta(i,l) * KroneckerDelta(j,k) ) ) ) ...
                               + ( S3(i) * KroneckerDelta(i,j) * Nhat(k) * Nhat(l) ) ... 
                               + ( S4(k) * KroneckerDelta(k,l) * Nhat(i) * Nhat(j) ) ...
                               + ( S5(i) * ( ( KroneckerDelta(i,k) * Nhat(j) * Nhat(l) ) ...
                                           + ( KroneckerDelta(i,l) * Nhat(j) * Nhat(k) ) ) ) ...
                               + ( S6(j) * ( ( KroneckerDelta(j,k) * Nhat(i) * Nhat(l) ) ...
                                           + ( KroneckerDelta(j,l) * Nhat(i) * Nhat(k) ) ) ) ...
                               + ( ArrayS7(i,j,k,l) * Nhat(i) * Nhat(j) * Nhat(k) * Nhat(l) ) ) ; 
                end ; 
            end ; 
        end ; 
    end ; 
       
else 
    
    %   sphere 
    %disp('Spherical inclusion...') ; 
    
    RhoSphere = a / sqrt(rSquared) ; 
    RhoSphereSquared = RhoSphere^2 ; 
    
    %   Ju & Sun 2001, eqn A.11 
    for i = 1:3 
        for j = 1:3 
            for k = 1:3 
                for l = 1:3 
                    ArrayGbar(i,j,k,l) = ( ( RhoSphere * RhoSphereSquared ) / ( 30 * OneMinusNu ) ) ...
                                           * ( ( ( 3 * RhoSphereSquared + 10 * Nu - 5 ) ...
                                            * KroneckerDelta(i,j) * KroneckerDelta(k,l) ) ...
                                            + ( ( 3 * RhoSphereSquared - 10 * Nu + 5 ) ...
                                            * ( ( KroneckerDelta(i,k) * KroneckerDelta(j,l) ) ...
                                              + ( KroneckerDelta(i,l) * KroneckerDelta(j,k) ) ) ) ...
                                            + ( ( 15 * ( 1 - RhoSphereSquared ) ) ...
                                               * KroneckerDelta(i,j) * Nhat(k) * Nhat(l) ) ...
                                            + ( ( 15 * ( 1 - Nu2 - RhoSphereSquared ) ) ...
                                                * KroneckerDelta(k,l) * Nhat(i) * Nhat(j) ) ...
                                            + ( ( 15 * ( Nu - RhoSphereSquared ) ) ... 
                                             * ( ( KroneckerDelta(i,k) * Nhat(j) * Nhat(l) ) ...
                                               + ( KroneckerDelta(i,l) * Nhat(j) * Nhat(k) ) ...
                                               + ( KroneckerDelta(j,k) * Nhat(i) * Nhat(l) ) ...
                                               + ( KroneckerDelta(j,l) * Nhat(i) * Nhat(k) ) ) ) ...
                                            + ( ( 15 * ( 7 * RhoSphereSquared - 5 ) ) ...
                                            * Nhat(i) * Nhat(j) * Nhat(k) * Nhat(l) ) ) ; 
                end ; 
            end; 
        end ;
    end ; 
        
end ;   

%   map the 3x3x3x3 tensor into 6x6 matrix, in Voigt order 
EshelbyG(1,1) = ArrayGbar(1,1,1,1) ; 
EshelbyG(2,2) = ArrayGbar(2,2,2,2) ; 
EshelbyG(3,3) = ArrayGbar(3,3,3,3) ; 
EshelbyG(4,4) = ArrayGbar(2,3,2,3) ; 
EshelbyG(5,5) = ArrayGbar(3,1,3,1) ; 
EshelbyG(6,6) = ArrayGbar(1,2,1,2) ; 
EshelbyG(1,2) = ArrayGbar(1,1,2,2) ; 
EshelbyG(1,3) = ArrayGbar(1,1,3,3) ; 
EshelbyG(2,1) = ArrayGbar(2,2,1,1) ; 
EshelbyG(3,1) = ArrayGbar(3,3,1,1) ; 
EshelbyG(3,2) = ArrayGbar(3,3,2,2) ; 
EshelbyG(2,3) = ArrayGbar(2,2,3,3) ; 