%   EshelbyDemo.m
%
%   Eshelby demo routine 
%
%   Calls various routines to calculate the elastic deformation
%   (i.e. stress and strain) in and around an ellipsoidal inclusion (or 
%   'heterogeneity') within an infinite matrix
%
%   Notes:
%       1. as a novice with MATLAB - you can probably guess! - I've
%       deliberately favoured an 'explicit' layout of the
%       various equations to make it easier to debug; 
%       I reckon these routines could be shorter and/or quicker; but as they are, I
%       can see how the statements map to the original published
%       equations...
%       2. semi-axis a // x-axis, b // y-axis, c // z-axis
%       3. ellipsoid centred at x = 0, y = 0, z = 0 
%       4. linear elastic isotropy in the inclusion and the matrix
%       
%   Based on the formulations originally derived by Eshelby:
%       Eshelby, 1957. Proc. Royal Society London. A241, p376. 
%       Eshelby, 1959. Proc. Royal Society London. A251, p561. 
%   and then clarified by Mura:
%       Mura, 1987. Micromechanics of Defects in Solids. 2nd ed. Dordrecht.
%   and Ju & Sun:
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

% clear all ; 
close all ; 

disp(' ') ; 
disp(['Starting EshelbyDemo at ', datestr(now), ',,,']) ; 
disp(' ') ; 

%   ellipsoid dimensions
a = 0.2 ; 
b = 1 ;
c = 1 ;
n = 100.0 ; 
lim = 4 ; 
[xe, ye, ze] = ellipsoid(0.0, 0.0, 0.0, a, b, c, n) ;

%   ellipsoid plot
figure ; 
set(gcf, 'PaperPositionMode', 'manual') ; 
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

surf(xe, ye, ze, 'FaceColor', 'red', 'EdgeColor','none') ; 

axis equal ; 
xlim([-lim lim]) ; 
ylim([-lim lim]) ; 
zlim([-lim lim]) ; 
set(gca, 'XTickMode', 'manual');
set(gca, 'YTickMode', 'manual');
set(gca, 'ZTickMode', 'manual');
set(gca,'XTick',[-2 0 2]) ; 
set(gca,'YTick',[-2 0 2]) ; 
set(gca,'ZTick',[-2 0 2]) ; 
grid on ;
camlight('left') ; 
lighting phong ; 
title('Ellipsoidal void geometry') ; 

legend('Ellipsoid', 'Location', 'southeast') ; 

print -r300 -dtiff 'Eshelby_ellipsoid.tiff'

%   1. input parameters 

%   ellipsoid elastic constants
Inc_E = 0 * 1e9 ; 
Inc_Nu = 0 ; 

%   matrix elastic constants 
Mat_E = 50.0 * 1e9 ; 
Mat_Nu = 0.25 ; 

%   calculate stiffness & compliance for inclusion 
Inc_Stiffness = calcIsotropicStiffness(Inc_E, Inc_Nu) ; 
% Inc_Compliance = inv(Inc_Stiffness) ; 

%   calculate stiffness & compliance for matrix 
Mat_Stiffness = calcIsotropicStiffness(Mat_E, Mat_Nu) ; 
% Mat_Compliance = inv(Mat_Stiffness) ; 

%   calculate S tensor 
Eshelby_S = calcEshelbyS(Mat_Nu, a, b, c) ; 

disp('Eshelby S tensor - ') ; 
disp(Eshelby_S) ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Run 1 - tensile stress normal to the crack 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   prescribed inclusion strain, column vector in Voigt order 
%   or 'transformation strain' 
%   epsilon_XX, epsilon_YY, epsilon_ZZ, gamma_YZ, gamma_XZ, gamma_XY
Inc_Strain = [ 0 ; ...
               0 ; ... 
               0 ; ... 
               0 ; ... 
               0 ; ... 
               0 ] ; 

%   prescribed matrix stress, column vector in Voigt order, in Pa 
%   sigma_XX, sigma_YY, sigma_ZZ, tau_YZ, tau_XZ, tau_XY
Mat_Stress = [ 1e6 ; ... 
               0 ; ... 
               0 ; ...
               0 ; ...
               0 ; ...
               0 ] ; 

%   calculate initial matrix strain from stress 
% Mat_Strain = Mat_Compliance * Mat_Stress ; 
Mat_Strain = Mat_Stiffness \ Mat_Stress ; 
% disp('Initial matrix strain - ') ; 
% disp(Mat_Strain) ; 
% disp('Initial matrix stress - ') ; 
% disp(Mat_Stress) ; 
% 
%   calculate total eigenstrain  
Eigen_Strain = calcEigenStrain(Mat_Stiffness, ... 
                               Inc_Stiffness, ... 
                               Eshelby_S, ...
                               Mat_Strain, ...
                               Inc_Strain) ;

%   calculate inclusion internal field (constant!)
Identity4 = eye(6, 6) ; 

%   Ju & Sun 1999, eqn 5 
Internal_Strain = Eshelby_S * Eigen_Strain ; 

%   Ju & Sun 1999, eqn 6
Internal_Stress = ( Mat_Stiffness * ( Eshelby_S - Identity4 ) ) * Eigen_Strain ; 

% disp('Total internal inclusion strain - ') ; 
% disp(Internal_Strain) ; 
% disp('Total internal inclusion stress - ') ; 
% disp(Internal_Stress) ; 
% 
%   calculate matrix total strain - the external field 
%   only calculate in the positive cartesian octant
%   with a >= b >= c and a // x-axis, b // y-axis, c // z-axis
Num_increments = 20 ; 
Max_x = lim ; 
Max_y = Max_x ; 
Max_z = Max_x ; 
Incr_x = Max_x / Num_increments ; 
Incr_y = Incr_x ;
Incr_z = Incr_x ; 

%   XY plane, z = 0 
xi = 0 ; 
Dim = 2 * Num_increments + 1 ; 
xcoord = zeros(Dim, Dim, Dim) ; 
ycoord = zeros(Dim, Dim, Dim) ; 
zcoord = zeros(Dim, Dim, Dim) ; 
rcoord = zeros(Dim, Dim, Dim) ; 
sigmaxx = zeros(Dim, Dim, Dim) ; 
sigmayy = zeros(Dim, Dim, Dim) ; 
sigmazz = zeros(Dim, Dim, Dim) ; 
sigmaxy = zeros(Dim, Dim, Dim) ; 
sigmaMean = zeros(Dim, Dim, Dim) ; 

for x = -Max_x:Incr_x:Max_x
    
    xi = xi + 1 ; 
    yi = 0 ; 
    
    for y = -Max_y:Incr_y:Max_y 
        
        yi = yi + 1 ; 
        zi = 0 ; 
        
        for z = -Max_z:Incr_z:Max_z 
        
            zi = zi + 1 ;         
        
            xcoord(xi,yi,zi) = x ; 
            ycoord(xi,yi,zi) = y ; 
            zcoord(xi,yi,zi) = z ; 

            %   calculate r, as the distance of each obs point from the origin 
            rcoord(xi, yi, zi) = sqrt(x^2 + y^2 + z^2) ; 
            
            %   only work out the external field for external points!
            if ( IsPointOutside(x, y, z, a, b, c) ) 
            
                External_Stress = calcExternalField(x, y, z, ... 
                                                    a, b, c, ... 
                                                    Mat_Nu, ...
                                                    Eigen_Strain, ...
                                                    Mat_Stiffness) ; 
                
                %   use Hooke's law to get external strain from stress 
%                 External_Strain = Mat_Compliance * External_Stress ; 
                External_Strain = Mat_Stiffness \ External_Stress ; 
                
                %   total field = initial matrix field + inclusion field  
                Total_Strain = External_Strain ; 
                Total_Stress = External_Stress ; 
                
                %   store stress 
                sigmaxx(xi,yi,zi) = Total_Stress(1) ; 
                sigmayy(xi,yi,zi) = Total_Stress(2) ; 
                sigmazz(xi,yi,zi) = Total_Stress(3) ; 
            
                sigmaxy(xi,yi,zi) = Total_Stress(6) ; 
            
                sigmaMean(xi,yi,zi) = ( sigmaxx(xi,yi,zi) + ...
                                        sigmayy(xi,yi,zi) + ... 
                                        sigmazz(xi,yi,zi) ) / 3 ; 
            
            else

                sigmaxx(xi,yi,zi) = NaN ; 
                sigmayy(xi,yi,zi) = NaN ; 
                sigmazz(xi,yi,zi) = NaN ; 
                sigmaxy(xi,yi,zi) = NaN ; 
            
                sigmaMean(xi,yi,zi) = NaN ; 
        
            end 
        
        end
        
    end
    
end

%   graph external field crack normal stress, in 3D  
figure ; 
set(gcf, 'PaperPositionMode', 'manual') ; 
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

surf(xe, ye, ze, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.4) ; 

ptens = patch(isosurface(xcoord,ycoord,zcoord,sigmaxx,+2.5e4)) ; 
set(ptens, 'FaceColor', [1,0.5,0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
p0 = patch(isosurface(xcoord,ycoord,zcoord,sigmaxx,0)) ; 
set(p0, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.4);
pcomp = patch(isosurface(xcoord,ycoord,zcoord,sigmaxx,-3e4)) ; 
set(pcomp, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.4);

axis equal ; 
xlim([-lim lim]) ; 
ylim([-lim lim]) ; 
zlim([-lim lim]) ; 
set(gca, 'XTickMode', 'manual');
set(gca, 'YTickMode', 'manual');
set(gca, 'ZTickMode', 'manual');
set(gca,'XTick',[-2 0 2]) ; 
set(gca,'YTick',[-2 0 2]) ; 
set(gca,'ZTick',[-2 0 2]) ; 
grid on ; 
daspect([1 1 1]) ; 
view(3); 
camlight ; 
lighting phong ; 
title('Stresses due to crack normal tension') ; 

[ ~, BLicons ] = legend('Ellipsoid', '\sigma_{XX}=+25 kPa', ...
                        '\sigma_{XX}=0 kPa', '\sigma_{XX}=-30 kPa', ...
                        'Location', 'southeast') ; 
PatchInLegend = findobj(BLicons, 'type', 'patch') ;
set(PatchInLegend, 'facea', 0.4) ; 

print -r300 -dtiff 'Eshelby_all_stress.tiff'

%   just tensile contours 
figure ; 
set(gcf, 'PaperPositionMode', 'manual') ; 
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

surf(xe, ye, ze, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.4) ; 

ptens2 = patch(isosurface(xcoord,ycoord,zcoord,sigmaxx,+4e4)) ; 
set(ptens2, 'FaceColor', [1,0.15,0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

ptens1 = patch(isosurface(xcoord,ycoord,zcoord,sigmaxx,+2.5e4)) ; 
set(ptens1, 'FaceColor', [1,0.5,0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

ptens3 = patch(isosurface(xcoord,ycoord,zcoord,sigmaxx,+2e4)) ; 
set(ptens3, 'FaceColor', [1,1,0], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

axis equal ; 
xlim([-lim lim]) ; 
ylim([-lim lim]) ; 
zlim([-lim lim]) ; 
set(gca, 'XTickMode', 'manual');
set(gca, 'YTickMode', 'manual');
set(gca, 'ZTickMode', 'manual');
set(gca,'XTick',[-2 0 2]) ; 
set(gca,'YTick',[-2 0 2]) ; 
set(gca,'ZTick',[-2 0 2]) ; 
grid on ; 
daspect([1 1 1]) ; 
view(3); 
camlight ; 
lighting phong ; 
title('Tensile stresses due to crack normal tension') ; 

[ ~, BLicons ] = legend('Ellipsoid', '\sigma_{XX}=+40 kPa', ...
                        '\sigma_{XX}=+25 kPa', '\sigma_{XX}=+20 kPa', ... 
                        'Location', 'southeast') ; 
PatchInLegend = findobj(BLicons, 'type', 'patch') ;
set(PatchInLegend, 'facea', 0.4) ; 

print -r300 -dtiff 'Eshelby_tensile_stress.tiff'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   calculate locus of max CNS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xcns = xcoord(find(~ycoord)) ; 
zcns = zcoord(find(~ycoord)) ; 
xzcns = sigmaxx(find(~ycoord)) ; 

xcns = reshape(xcns, [Dim, Dim]) ; 
zcns = reshape(zcns, [Dim, Dim]) ; 
xzcns = reshape(xzcns, [Dim, Dim]) ; 

figure ; 
set(gcf, 'PaperPositionMode', 'manual') ; 
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

hold on ; 
contourf(xcns, zcns, xzcns, 15) ; 
contour(xcns, zcns, xzcns, [0 0]) ; 
hold off ; 

axis equal ; 
xlim([0 lim]) ; 
ylim([0 lim]) ; 
grid on ; 
box on ; 
title('XZ CNS') ; 
colorbar ; 
% caxis([-2e6 2e6]) ; 

xcns = xcoord(find(~zcoord)) ; 
ycns = ycoord(find(~zcoord)) ; 
xycns = sigmaxx(find(~zcoord)) ; 

xcns = reshape(xcns, [Dim, Dim]) ; 
ycns = reshape(ycns, [Dim, Dim]) ; 
xycns = reshape(xycns, [Dim, Dim]) ; 

figure ; 
set(gcf, 'PaperPositionMode', 'manual') ; 
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

hold on ; 
contourf(xcns, ycns, xycns, 15) ; 
contour(xcns, ycns, xycns, [0 0]) ; 
hold off ; 

axis equal ; 
xlim([0 lim]) ; 
ylim([0 lim]) ; 
grid on ; 
box on ; 
title('XY CNS') ; 
colorbar ; 
% caxis([-2e6 2e6]) ; 

%   clean-up & quit 
disp(' ') ; 
disp(['...finished EshelbyDemo at ', datestr(now), '.']) ; 
disp(' ') ; 
