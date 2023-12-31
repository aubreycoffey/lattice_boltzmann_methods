%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lecture on Lattice Boltzmann methods                                %%%
%%% TU M�nchen summer term 2021                                         %%%
%%%                                                                     %%%
%%% Prof. Dr. Barbara Wohlmuth                                          %%%
%%% M.Sc. Gladys GUtierrez                                              %%%
%%% M.Sc. Markus Muhr                                                   %%%
%%%                                                                     %%%
%%% Programming Tutorial 1 - HPP Implementation                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Channel geometry
Nx = 200;
Ny = 100;
Nt = 1000;
plot_interval = 100;

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

% Logical array indicating obstacle cells
% square
obstacle = [cyl_y,cyl_x];
for i=cyl_x-cyl_r: cyl_x+cyl_r
    for j=cyl_y-cyl_r: cyl_y+cyl_r
        dist=(((cyl_x-i)^2)+((cyl_y-j)^2))^.5;
        if dist<=cyl_r
            obstacle=[obstacle;[j,i]];
        end
    end
end
%use eq of a circle to define obstacle, use for loops
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise b)                                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify boundary cells
bbcells = zeros(Ny,Nx);
% For which cells are the bounce back rules applied? When exercise c) is
% done, don't forget to take into account the obstacle created 
bbcells(Ny,:) = 1;          
bbcells(1,:) = 1;             

%set obstacle to 0
%bbcells(obstacle)=0;
for k = 1:size(obstacle,1)
    a=obstacle(k,:)
    bbcells(obstacle(k,:))=1;
end
%top and bottom
bbcells = find(bbcells == 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%random number(0,1) vec from ny to 1, check if its smaller than velocity profile


%%% Graphical sketch of the numbering in the HPP stencil
%%%
%%%    3     2
%%%     \   /
%%%   4  - o - 1
%%%     /   \
%%%    5     6

% Discrete velocities
% cx  = 1/sqrt(2) * [1, 1, -1, -1,-1,1];
% cy  = 1/sqrt(2) * [0, 1, 1, 0,-1,-1];
%above is first attempt then remembered the lattice is equilateral
%triangles so we'll have to use 30-60-90 triangles

%for the diagonal pieces if the horizontal is 1 the vertical is 2 and the
%norm is sqrt(3)
% cx  = 1/sqrt(3) * [sqrt(3), 1, -1, -sqrt(3),-1,1];
% cy  = 1/sqrt(3) * [0, 2, 2, 0,-2,-2];
% opposing directions (necessary for boundary treatment)
opp = [ 4,5,6, 1, 2,3 ];

% Initialization of population arrays
fIn  = zeros(6,Ny,Nx);
fOut = fIn;

% Random (0,1) Initial data
for i = 1:6
    fIn(i,:,:) = randi([0 1],[Ny,Nx]); %.. .. ..
end

% Main time loop
for t = 1:Nt

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Exercise b)                                                     %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set boundary condition at Inlet, probability is given by the
    %normalized Poiseuille profile, use logical operator to compare
    %random generated values between 0 and 1 and the velocity profile
    %values
    for i = [1,6]%
        %fIn(i,2:Ny-1,1) = round(vel_profile([2:Ny-1])) %gave us i
        fIn(i,2:Ny-1,1) = rand(1,Ny-2,1)<= vel_profile([2:Ny-1])
    end
    
    % Set boundary condition at Outlet
    for i = [2,3]%
        fIn(i,:,Nx) = 0;%gave us i
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Macroscopic quantities
    rho = sum(fIn,1); %fIn(i,:,:) %.. .. .. fin summation
    ux  = reshape ( (cx * reshape(fIn,6,Ny*Nx)), 1,Ny,Nx) ./rho;
    ux(isnan(ux)) = 0; %.. .. ..
    uy  = reshape ( (cy * reshape(fIn,6,Ny*Nx)), 1,Ny,Nx) ./rho; %.. .. ..
    uy(isnan(uy)) = 0; %.. .. ..

    % Collision step
    for i = 1:6        
        % vector that indicates counter clockwise HPP numbering
        rf = [ 1, 2, 3, 4,5,6, 1, 2, 3, 4,5,6 ];
        
        % Collision function (see exercise)
        %Delta = -fIn(rf(i),:,:) .* (1 - fIn(rf(i+1),:,:)) .*.. .. ..
        Delta = (-fIn(rf(i),:,:) .* (1 - fIn(rf(i+1),:,:)) .* fIn(rf(i+2),:,:).* (1 - fIn(rf(i+3),:,:)))+((1 - fIn(rf(i),:,:)) .*fIn(rf(i+1),:,:) .*(1 - fIn(rf(i+2),:,:)).*fIn(rf(i+3),:,:));
        %?i(n) = ? ni(1 ? ni+1)ni+2(1 ? ni+3)+ (1 ? ni)ni+1(1 ? ni+2)ni+3

        rv=randi([0 1],1)
        Delta = (fIn(rf(i+1),:,:) .* fIn(rf(i+3),:,:) .* fIn(rf(i+5),:,:).* (1 - fIn(rf(i),:,:)).* (1 - fIn(rf(i+2),:,:)).* (1 - fIn(rf(i+4),:,:)))...
        -(fIn(rf(i),:,:) .* fIn(rf(i+2),:,:) .* fIn(rf(i+4),:,:).* (1 - fIn(rf(i+1),:,:)).* (1 - fIn(rf(i+3),:,:)).* (1 - fIn(rf(i+5),:,:)))...
        +(rv*fIn(rf(i+1),:,:).* fIn(rf(i+4),:,:).* (1 - fIn(rf(i),:,:)).* (1 - fIn(rf(i+2),:,:)).* (1 - fIn(rf(i+3),:,:)).* (1 - fIn(rf(i+5),:,:)))...
        +((1-rv)*fIn(rf(i+2),:,:).* fIn(rf(i+5),:,:).* (1 - fIn(rf(i),:,:)).* (1 - fIn(rf(i+1),:,:)).* (1 - fIn(rf(i+3),:,:)).* (1 - fIn(rf(i+4),:,:)))...
        -(fIn(rf(i),:,:) .* fIn(rf(i+3),:,:).* (1 - fIn(rf(i+1),:,:)).* (1 - fIn(rf(i+2),:,:)).* (1 - fIn(rf(i+4),:,:)).* (1 - fIn(rf(i+5),:,:)));
      
        
        % Population update
        fOut(i,:,:) = fIn(i,:,:)+Delta(i);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Exercise b)                                                     %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  use opp vector created previously. Applying bounce back means having a
%  zero velocity at the walls, we achieve this by setting the sum of c_i*f_i
%  equal to zero: this is done when probailities of opposite directions velocities are the same 
    for i = 1:6
         fOut(i,bbcells) = fIn(opp(i),bbcells);%
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Streaming step
    for i = 1:6
        fIn(i,:,:) = circshift(fOut(i,:,:),[0,-sqrt(2)*cy(i),sqrt(2)*cx(i)]);
        %[0,-cy(i),cx(i)] might need a sqrt2 infront of cy cx
        %circshift(fOut(i,:,:),.. .. .., .. .. ..);
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
                    %if n==0 there are no points in the average because
                    %were in the obstacle, im guessing? (since n==0 theres 
                    %nothing being averaged, ie no neighbors) then we want
                    %ux,uy== 1 or 0, or it will accomplish that by letting
                    %it sum
                    ux(1,i,j) = sum(sum(UX(range_y,range_x) .* ~obstacle(range_y,range_x),1),2);
                    uy(1,i,j) = sum(sum(UY(range_y,range_x) .* ~obstacle(range_y,range_x),1),2);%
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
        U = (UX.^2+UY.^2 ).^(.5)% .. .. ..
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






