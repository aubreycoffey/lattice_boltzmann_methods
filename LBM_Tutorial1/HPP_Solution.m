%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lecture on Lattice Boltzmann methods                                %%%
%%% TU München summer term 2017                                         %%%
%%%                                                                     %%%
%%% Prof. Dr. Barbara Wohlmuth                                          %%%
%%% M.Sc. Gladys Gutierrez                                              %%%
%%% M.Sc. Markus Muhr                                                   %%%
%%%                                                                     %%%
%%% Programming Tutorial 1 - HPP Implementation                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Channel geometry
Nx = 200;
Ny = 100;
Nt = 1000;
plot_interval = 500;

% Cell coordinates
[y,x] = meshgrid(1:Ny,1:Nx); 
y     = flipud(y');                        
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


%%% Graphical sketch of the numbering in the HPP stencil
%%%
%%%    2     1
%%%     \   /
%%%       o
%%%     /   \
%%%    3     4

% Discrete velocities
cx  = 1/sqrt(2) * [ 1, -1, -1,  1 ];
cy  = 1/sqrt(2) * [ 1,  1, -1, -1 ];

% opposing directions (necessary for boundary treatment)
opp = [ 3, 4, 1, 2 ];

% Initialization of population arrays
fIn  = zeros(4,Ny,Nx);
fOut = fIn;

% Initial data
for i = 1:4
    fIn(i,:,:) = randi(2,Ny,Nx) - 1;       % Exercise a)
    fIn(i,:,:) = 0;                        % Exercise b) following
end

% Main time loop
for t = 1:Nt

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Exercise b)                                                     %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set boundary condition at Inlet
    for i = [1,4]
        fIn(i,[2:Ny-1],1) = rand(1,Ny-2,1) <= vel_profile([2:Ny-1]);
    end
    for i = [2,3]
        fIn(i,[2:Ny-1],1) = 0;
    end
    
    % Set boundary condition at Outlet
    for i = [2,3]
        fIn(i,[2:Ny-1],Nx) = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Macroscopic quantities
    rho = sum(fIn,1);
    ux  = reshape ( (cx * reshape(fIn,4,Ny*Nx)), 1,Ny,Nx) ./rho;
    ux(isnan(ux)) = 0;
    uy  = reshape ( (cy * reshape(fIn,4,Ny*Nx)), 1,Ny,Nx) ./rho; 
    uy(isnan(uy)) = 0;

    % Collision step
    for i = 1:4        
        % vector that indicates counter clockwise HPP numbering
        rf = [ 1, 2, 3, 4, 1, 2, 3, 4 ];
        
        % Collision function (see exercise)
        Delta = -fIn(rf(i),:,:) .* (1 - fIn(rf(i+1),:,:)) .* fIn(rf(i+2),:,:) .* (1 - fIn(rf(i+3),:,:))...
                +(1 - fIn(rf(i),:,:)) .* fIn(rf(i+1),:,:) .* (1 - fIn(rf(i+2),:,:)) .* fIn(rf(i+3),:,:);
        
        % Population update
        fOut(i,:,:) = fIn(i,:,:) + Delta;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Exercise b)                                                     %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:4
         fOut(i,bbcells) = fIn(opp(i),bbcells);
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Streaming step
    for i = 1:4
        fIn(i,:,:) = circshift(fOut(i,:,:), [0,-sqrt(2)*cy(i),sqrt(2)*cx(i)]);
    end
    
    
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
        M = zeros(2*Ny,2*Nx);
        M([1:2:2*Ny-1],[1:2:2*Nx-1]) = fIn(2,:,:);
        M([2:2:2*Ny],[1:2:2*Nx-1]) = fIn(3,:,:);
        M([1:2:2*Ny-1],[2:2:2*Nx]) = fIn(1,:,:);
        M([2:2:2*Ny],[2:2:2*Nx]) = fIn(4,:,:);
        M([1,2,2*Ny-1,2*Ny],:) = 0;

        figure(1);
        spy(M);
        tt = title('HPP','Interpreter','LaTex');
        set(tt,'FontSize',20);
        ax = handle(gca);
        ax.XTick = 0:1:2*Nx;
        ax.YTick = 0:1:2*Ny;
        ax.YTickLabel = [];
        ax.XTickLabel = [];
        ax.GridLineStyle = '-';
        ax.XColor = [0.3,0.3,0.3];
        ax.YColor = [0.3,0.3,0.3];
        xlabel('');
        ylim([1, 2*Ny]);
        xlim([1, 2*Nx]);
        grid on;
        pause(0.1);
        drawnow;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Exercise b)                                                 %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
end






