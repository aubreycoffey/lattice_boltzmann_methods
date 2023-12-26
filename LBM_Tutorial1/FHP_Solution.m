%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lecture on Lattice Boltzmann methods                                %%%
%%% TU München summer term 2017                                         %%%
%%%                                                                     %%%
%%% M.Sc. Gladys Gutierrez                                              %%%
%%% Dr. Laura Scarabosio                                                %%%
%%% M.Sc. Markus Muhr                                                   %%%
%%%                                                                     %%%
%%% Programming Tutorial 1 - FHP-I Implementation                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Channel geometry
Nx = 200;
Ny = 100;
Nt = 1000;
plot_interval = 500;

% Cell coordinates
[y,x] = meshgrid(1:Ny,1:Nx); 
y     = flipud(y');    
y(:,[1:2:Nx]) = y(:,[1:2:Nx]) + 0.5;    % <-------------------- This is new
x     = x';

% Inflow profile (Poiseuille normed to max velocity of 1)
vel_profile = @(y) 1 / (Ny/2-1/2)^2 * (y-1).*(Ny - y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Exercise c)                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obstacle configuration
cyl_x    = Nx/5 + 1;                      % x-position of cylinder midpoint
cyl_y    = Ny/2;                          % y-position of cylinder midpoint
cyl_r    = Ny/10 + 1;                     % radius of cylinder
obstacle = (x - cyl_x).^2 + (y - cyl_y).^2 <= cyl_r.^2 - 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Exercise b)                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify boundary cells
bbcells = zeros(Ny,Nx);
bbcells(obstacle) = 1;                              % Exercise c) following
bbcells([1,Ny],:) = 1;             
bbcells = find(bbcells == 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Graphical sketch of the numbering in the FHP-I stencil
%%%
%%%   2     1
%%%    \   /
%%% 3 -- 0 -- 6
%%%    /   \
%%%   4     5

% Discrete velocities    <-------------------- This is new (new velocities)
cx  = [ cos(pi/3), cos(2*pi/3), cos(pi),  cos(4*pi/3), cos(5*pi/3), cos(2*pi) ];
cy  = [ sin(pi/3), sin(2*pi/3), sin(pi),  sin(4*pi/3), sin(5*pi/3), sin(2*pi) ];

% opposing directions (necessary for boundary treatment)
opp = [ 4, 5, 6, 1, 2, 3 ];             % <-------------------- This is new

% Initialization of population arrays
fIn  = zeros(6,Ny,Nx);  % <-------------------- This is new (array size is now 6)
fOut = fIn;

% Initial data
for i = 1:6        % <-------------------- This is new (loop size is now 6)
    fIn(i,:,:) = randi(2,Ny,Nx) - 1;       % Exercise a)
    fIn(i,:,:) = 0;                        % Exercise b) following
end

% Main time loop
for t = 1:Nt

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Exercise b)                                                     %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set boundary condition at Inlet
    for i = [1,5,6]                     % <-------------------- This is new
        fIn(i,[2:Ny-1],1) = rand(1,Ny-2,1) <= vel_profile([2:Ny-1]);
    end
    for i = [2,3,4]                     % <-------------------- This is new
        fIn(i,[2:Ny-1],1) = 0;
    end
    
    % Set boundary condition at Outlet
    for i = [2,3,4]                     % <-------------------- This is new
        fIn(i,[2:Ny-1],Nx) = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Macroscopic quantities
    rho = sum(fIn,1);
    ux  = reshape ( (cx * reshape(fIn,6,Ny*Nx)), 1,Ny,Nx) ./rho; % <-- (4 --> 6)
    ux(isnan(ux)) = 0;
    uy  = reshape ( (cy * reshape(fIn,6,Ny*Nx)), 1,Ny,Nx) ./rho; % <-- (4 --> 6)
    uy(isnan(uy)) = 0;

    % Collision step
    
    % Random variable for two-particle head on collision
        xi = randi(2,1,Ny,Nx) - 1;      % <-------------------- This is new
        
    for i = 1:6    % <-------------------- This is new (loop size is now 6)        
        % vector that indicates counter clockwise FHP numbering
        rf = [ 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6 ]; % <------- This is new
        
        % Collision function (see exercise) % <----- New collision function
        Delta = fIn(rf(i+1),:,:) .* fIn(rf(i+3),:,:) .* fIn(rf(i+5),:,:) .* (1 - fIn(rf(i),:,:)) .* (1 - fIn(rf(i+2),:,:)) .*(1 - fIn(rf(i+4),:,:)) +...
                -fIn(rf(i),:,:) .* fIn(rf(i+2),:,:) .* fIn(rf(i+4),:,:) .* (1 - fIn(rf(i+1),:,:)) .* (1 - fIn(rf(i+3),:,:)) .* (1 - fIn(rf(i+5),:,:)) +...
                +xi .* fIn(rf(i+1),:,:) .* fIn(rf(i+4),:,:) .* (1 - fIn(rf(i),:,:)) .* (1 - fIn(rf(i+2),:,:)) .* (1 - fIn(rf(i+3),:,:)) .* (1 - fIn(rf(i+5),:,:)) +...
                +(1 - xi) .* fIn(rf(i+2),:,:) .* fIn(rf(i+5),:,:) .* (1 - fIn(rf(i),:,:)) .* (1 - fIn(rf(i+1),:,:)) .* (1 - fIn(rf(i+3),:,:)) .* (1 - fIn(rf(i+4),:,:)) + ...
                - fIn(rf(i),:,:) .* fIn(rf(i+3),:,:) .* (1 - fIn(rf(i+1),:,:)) .* (1 - fIn(rf(i+2),:,:)) .* (1 - fIn(rf(i+4),:,:)) .* (1 - fIn(rf(i+5),:,:));
        
        % Population update
        fOut(i,:,:) = fIn(i,:,:) + Delta;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Exercise b)                                                     %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:6    % <-------------------- This is new (loop size is now 6)
         fOut(i,bbcells) = fIn(opp(i),bbcells);
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Streaming step                    % <-------------------- This is new
    fIn(1,:,:) = circshift(fOut(1,:,:), [0,-1,1]);
    fIn(2,:,:) = circshift(fOut(2,:,:), [0,-1,-1]);
    fIn(3,:,:) = circshift(fOut(3,:,:), [0,0,-2]);
    fIn(4,:,:) = circshift(fOut(4,:,:), [0,1,-1]);
    fIn(5,:,:) = circshift(fOut(5,:,:), [0,1,1]);
    fIn(6,:,:) = circshift(fOut(6,:,:), [0,0,2]);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Exercise d) Coarse graining                                     %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if t == Nt
        fw = 4;
        UX = zeros(Ny,Nx);
        UY = zeros(Ny,Nx);
        UX(:,:) = ux(1,:,:);
        UY(:,:) = uy(1,:,:);
        for i = 1:Ny
            for j = 1:Nx
                
                % Cut off the filter at the channel boundaries since you
                % dont want to filter with non existing cells
                range_y = [max(1,i-fw):min(Ny,i+fw)];
                range_x = [max(1,j-fw):min(Nx,j+fw)];
                
                % Also cut off the filter at the obstacle since there is no
                % flow inside and hence no cells to filter with
                n  = length(range_x)*length(range_y) - sum(sum(obstacle(range_y,range_x),1),2);
                if n==0
                    ux(1,i,j) = 0;
                    uy(1,i,j) = 0;
                else
                    ux(1,i,j) = 1/n*sum(sum(UX(range_y,range_x) .* ~obstacle(range_y,range_x),1),2);
                    uy(1,i,j) = 1/n*sum(sum(UY(range_y,range_x) .* ~obstacle(range_y,range_x),1),2);
                end
            end 
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % Graphical output
    if (mod(t,plot_interval) == 0) || (mod(t,plot_interval) == 1) || t == Nt
        figure(2)
        UX = zeros(Ny,Nx);
        UX(:,:) = ux(1,:,:);
        UX(bbcells) = 0;           
        UY = zeros(Ny,Nx);
        UY(:,:) = uy(1,:,:);
        UY(bbcells) = 0;           
        U = sqrt(UX.^2 + UY.^2);
        pcolor(flipud(U));
        shading flat;
        axis equal off; 
        tt = title('Velocity magnitude $|u|$','Interpreter','LaTex');
        set(tt,'FontSize',20);
        colorbar;
        colormap jet;
        drawnow;

        % VTK Output only at last time step to increase speed
        if t==Nt
            mkdir('VTKOutput');
            str = strcat('VTKOutput/velocity',num2str(t),'.vtk');
            vtkwrite(str,'structured_grid',x,y,zeros(size(x)),'vectors','velocity',UX,UY,zeros(size(UX)));
        end

    end
    
end
    







