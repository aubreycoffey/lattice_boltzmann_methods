%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lecture on Lattice Boltzmann methods                                %%%
%%% TU München summer term 2017                                         %%%
%%%                                                                     %%%
%%% Prof. Dr. Barbara Wohlmuth                                          %%%
%%% M.Sc. Gladys Gutierrez                                              %%%
%%% M.Sc. Markus Muhr                                                   %%%
%%%                                                                     %%%
%%% Programming Tutorial 2 - Carman-Vortex-Street D2Q9 LBM Simulation   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lattice Boltzmann sample in Matlab
% Copyright (C) 2006-2008 Jonas Latt
% Address: EPFL, 1015 Lausanne, Switzerland
% E-mail: jonas@lbmethod.org
% Get the most recent version of this file on LBMethod.org:
%   http://www.lbmethod.org/_media/numerics:cylinder.m
%
% Original implementaion of Zou/He boundary condition by
% Adriano Sciacovelli (see example "cavity.m")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public 
% License along with this program; if not, write to the Free 
% Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
% Boston, MA  02110-1301, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% DISCRETIZATION PARAMETERS (CHANNEL GEOMETRY)
Nx      = 400;                  % number of cells in x-direction
Ny      = 100;                  % number of cells in y-direction
Nt      = 5000;                 % number of time steps

[y,x]   = meshgrid(1:Ny,1:Nx);  % (x,y)-coordinates of each cell in cell units (this is helpful)  
y = flipud(y');                 % In order to match the x-y-orientation
x = x';

plot_interval = 150;            % How many timesteps for the next plot update


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Exercise a)                                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (INITIAL) OBSTACLE CONFIGURATION
cyl_x = Nx/5 + 1;               % x-position of cylinder midpoint
cyl_y = Ny/2 + 3;               % y-position of cylinder midpoint (exact y-symmetry is avoided)
cyl_r = Ny/10 + 1;              % radius of cylinder
obstacle = (x - cyl_x).^2 + (y - cyl_y).^2 <= cyl_r.^2 - 1; % Logical matrix indicating the obstacle cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Exercise b)                                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOUNCE BACK BOUNDARY CONDITION AT CHANNEL WALLS AND OBSTACLE   
%BBcells = .. .. ..

BBcells = zeros(Ny,Nx);
BBcells(obstacle) = 1;                              % Exercise c) following
BBcells([1,Ny],:) = 1;             
BBcells = find(BBcells == 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Exercise c)                                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLOW PARAMETERS
Re      = 100;                      % Reynolds number (Main flow parameter)0000
uMax    = 0.1;                      % Maximal velocity at channel inflow
%re=u*l/v, where l is diameter of cylinder
%v=u*l/re=u*(cyl_r*2)/re
nu      = (uMax*(cyl_r*2))/Re;  % kinematic viscosity (is chosen such that the given Raynolds number is achieved) 
%3? +1/2= ?, ?=1/omega
omega   = 1/((3*nu)+.5);         % relaxation parameter 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                          


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Exercise d)                                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   6   2   5
%%%     \ | /
%%%   3 - 0 - 1
%%%     / | \
%%%   7   4   8
% D2Q9 LATTICE CONSTANTS
w   = [4/9, 1/9, 1/9 , 1/9 , 1/9 , 1/36, 1/36 , 1/36 , 1/36 ];   % weights
cx  = [  0,   1, 0 , -1 , 0 ,    1, -1 , -1 , 1 ];   % x-coord of lattice velocities
cy  = [  0,   0, 1 , 0 , -1 ,    1, 1 , -1, -1];   % y-coord of lattice velocities
opp = [  1,   4, 5 , 2 , 3 ,  8 , 9 , 6 , 7 ];   % number of lattice velocity in opposite direction each
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                      

inlet  = 1;         % x-position of inlet (x-nr. of first cell)
outlet = Nx;        % x-position of outlet (x-nr. of last cell)

% Setting up iniital populations according to Poiseuille equilibrium:
%    Input distributions fIn(lattice velocity, cell x-coord., cell y-coord.)
%    Output distributions fOut (post collision see collision step later)
%    Equilibrium distribution fEq
fIn  = zeros(9,Ny,Nx);
fOut = fIn;
fEq  = fIn;



% MACROSCOPIC INITIAL CONDITION: Poiseuille profile at equilibrium
L      = Ny - 2;         % Number of fluid cells in y-direction (the bounce back cells at the channel walls do not count here)
y_phys = y - 1.5;        % Physical y-coordinates of fluid cells in order to evaluate the symmetric Poiseuille profile                
ux     = 4 * uMax / (L*L) * (y_phys.*L-y_phys.*y_phys); % x-veloxity according to Poiseuille   
uy     = zeros(Ny,Nx);                                  % y-velocity = 0 for all cells                                               
rho    = 1;              % Initial value for fluid density is 1 for every cell


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Exercise e)                                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MICROSCOPIC INITIAL CONDITION:
for i = 1:9               
    cu = 3 * ( cx(i)*ux + cy(i)*uy );
    un = sqrt(ux.^2 + uy.^2);
    fIn(i,:,:) = rho .* w(i) .* (1+cu+(.5*(cu.^2))+(-1.5*(un.^2)));    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% MAIN LOOP (TIME STEPS)
for t_step = 1:Nt
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Exercise f)                                                 %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE MACROSCOPIC QUANTITIES FOR EACH CELL
    rho = sum(fIn,1);
    ux  = reshape ( (cx * reshape(fIn,9,Ny*Nx)), 1,Ny,Nx) ./rho; %
    ux(isnan(ux)) = 0;
    uy  = reshape ( (cy * reshape(fIn,9,Ny*Nx)), 1,Ny,Nx) ./rho; % 
    uy(isnan(uy)) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % MACROSCOPIC (DIRICHLET) BOUNDARY CONDITIONS FOR THE FLOW
    
        % Inlet: Poiseuille profile for u (rho is set to a matching value)
            y_phys = [2:Ny-1] - 1.5;     % Again centering the inlet cells in y-direction  
            ux(:,[2:Ny-1],inlet) = 4 * uMax / (L*L) * (y_phys.*L-y_phys.*y_phys);
            uy(:,[2:Ny-1],inlet) = 0;        
            % rho is chosen in such a way that the velocity u is achieved
            % when using the microscopic Zhu/He BCs as follows
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Exercise g)                                                 %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            rho(:,[2:Ny-1],inlet) = ((1-ux(:,[2:Ny-1],inlet)).^(-1)).*(fIn(1,[2:Ny-1],inlet)+fIn(3,[2:Ny-1],inlet)+fIn(5,[2:Ny-1],inlet)+(2*(fIn(4,[2:Ny-1],inlet)+fIn(7,[2:Ny-1],inlet)+fIn(8,[2:Ny-1],inlet))));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % Outlet: Constant pressure (P = rho*R*T) hence we can also use constant density rho = 1 (u is set to a matching value)   
            rho(:,[2:Ny-1],outlet) = 1;        
            % Here rho is given, hence u has to be chosen in such a way,
            % that the microscopic Zhu/He conditions add up to the desired
            % value of rho
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Exercise g)                                                 %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
            ux(:,[2:Ny-1],outlet) = -1+((1./rho(:,[2:Ny-1],outlet)).*(fIn(1,[2:Ny-1],outlet)+fIn(3,[2:Ny-1],outlet)+fIn(5,[2:Ny-1],outlet)+(2*(fIn(2,[2:Ny-1],outlet)+fIn(6,[2:Ny-1],outlet)+fIn(9,[2:Ny-1],outlet)))));            
            uy(:,[2:Ny-1],outlet) = 0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   
            
    % MICROSCOPIC BOUNDARY CONDITIONS (Zou/He BC) FOR THE FLOW
    
        % Inlet: 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Exercise g)                                                 %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fIn(2,[2:Ny-1],inlet) = fIn(4,[2:Ny-1],inlet) + ((2/3)*rho(:,[2:Ny-1],inlet).*ux(:,[2:Ny-1],inlet));
            fIn(6,[2:Ny-1],inlet) = fIn(8,[2:Ny-1],inlet) + ((-.5)*(fIn(3,[2:Ny-1],inlet)-fIn(5,[2:Ny-1],inlet)))+ ((.5)*rho(:,[2:Ny-1],inlet).*uy(:,[2:Ny-1],inlet))+ ((1/6)*rho(:,[2:Ny-1],inlet).*ux(:,[2:Ny-1],inlet));
            fIn(9,[2:Ny-1],inlet) = fIn(7,[2:Ny-1],inlet) + ((.5)*(fIn(3,[2:Ny-1],inlet)-fIn(5,[2:Ny-1],inlet)))+ ((-.5)*rho(:,[2:Ny-1],inlet).*uy(:,[2:Ny-1],inlet))+ ((1/6)*rho(:,[2:Ny-1],inlet).*ux(:,[2:Ny-1],inlet));

        % Outlet:
            fIn(4,[2:Ny-1],outlet) = fIn(2,[2:Ny-1],outlet) + ((-2/3)*rho(:,[2:Ny-1],outlet).*ux(:,[2:Ny-1],outlet));
            fIn(8,[2:Ny-1],outlet) = fIn(6,[2:Ny-1],outlet) + ((-.5)*(fIn(5,[2:Ny-1],outlet)-fIn(3,[2:Ny-1],outlet)))+ ((-.5)*rho(:,[2:Ny-1],outlet).*uy(:,[2:Ny-1],outlet))+ ((-1/6)*rho(:,[2:Ny-1],outlet).*ux(:,[2:Ny-1],outlet));
            fIn(7,[2:Ny-1],outlet) = fIn(9,[2:Ny-1],outlet) + ((.5)*(fIn(5,[2:Ny-1],outlet)-fIn(3,[2:Ny-1],outlet)))+ ((.5)*rho(:,[2:Ny-1],outlet).*uy(:,[2:Ny-1],outlet))+ ((-1/6)*rho(:,[2:Ny-1],outlet).*ux(:,[2:Ny-1],outlet));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                                              
    % COLLISION STEP
    for i = 1:9
       cu = 3 * ( cx(i)*ux + cy(i)*uy );
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%% Exercise h)                                                 %%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Local equilibrium distributions
       un = sqrt(ux.^2 + uy.^2);
       fEq(i,:,:)  = rho(:,:,:) .* w(i) .* (1+cu+(.5*(cu.^2))+(-1.5*(un.^2)));    
       
       % Post collision distributions
       fOut(i,:,:) = fIn(i,:,:)-(omega*(fIn(i,:,:)-fEq(i,:,:)));
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    
    % OBSTACLE AND CHANNEL BOUNDARIES (BOUNCE-BACK)
    for i = 1:9
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%% Exercise i)                                                 %%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         fOut(i,BBcells) = fIn(opp(i),BBcells);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % STREAMING STEP
    for i = 1:9
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%% Exercise j)                                                 %%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       fIn(i,:,:) = circshift(fOut(i,:,:), [0,-cy(i),cx(i)]);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    % VISUALIZATION
    if (mod(t_step,plot_interval) == 1) 
        UX = zeros(Ny,Nx);
        UX(:,:) = ux(1,:,:);
        UX(BBcells) = NaN;           
        UY = zeros(Ny,Nx);
        UY(:,:) = uy(1,:,:);
        UY(BBcells) = NaN;   
        U = sqrt(UX.^2 + UY.^2);       
        pcolor(flipud(U));
        shading flat;
        axis equal off; 
        t = title('Velocity magnitude $|u|$','Interpreter','LaTex');
        set(t,'FontSize',20);
        colorbar;
        colormap jet;
        drawnow;          
    end
    
    % VTK Output only at last time step to increase speed
    if t_step == Nt
        UX(BBcells) = 0;
        UY(BBcells) = 0;
        str = strcat('VTKOutput/velocity',num2str(t_step),'.vtk');
        vtkwrite(str,'structured_grid',x,y,zeros(size(x)),'vectors','velocity',UX,UY,zeros(size(UX)));
    end
      
end