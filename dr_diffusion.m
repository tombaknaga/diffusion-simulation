%% My Fluid Sim
% Super simple fluid sim based on the diffusion equation
% du/dt = visc*(d2u/dx2 + d2u/dy2)
% u(x,y,t) is the concentration

% Spatial discretization using Central Difference

% Time discretization using the Explicit Euler method
% u_n+1 = u_n + h*f(t_n,y_n)

format compact
format short
clear all


%% Setting up the domain
lx = 0.02; % [m] length of the domain
ly = 0.02; % [m] height of the domain

% spatial discretization (mesh creation)
meshsize = 0.0001; % [m] mesh size
meshResX = lx/meshsize;
meshResY = ly/meshsize;

x = linspace(0,lx,meshResX);
y = linspace(0,ly,meshResY);

[X, Y] = meshgrid(x,y);

%% General settings
mu = 0.001; % [Pa.s] viscosity constant
rho = 1000; % [kg/m3] density
nu = mu/rho; % kinematic viscosity -> the constant in the diffusion eq.


dt = 0.1 * meshsize^2 / nu;
maxt = 5; % [s] maximum simulation time (real time)
t = linspace(0,maxt,maxt/dt);

%% Innitial condition
% To initialize the flow, we use the gaussian distribution centered at
% (lx/2, ly/2)
maxConMult = 1;

x0 = lx / 2; % Center x
y0 = ly / 2; % Center y
sigma = 0.001; % Standard deviation (spread of Gaussian)
u = maxConMult*exp(-((X - x0).^2 + (Y - y0).^2) / (2 * sigma^2)); % Initial concentration field

% save config data
save('initConfigData.mat',"X","Y","meshsize","mu","rho","maxt","dt","maxConMult","sigma");

if ~exist('dat', 'dir')
    mkdir('dat');
end

%% Into the time-loop!

% Time stepping loop
for it = 1:length(t)
    fprintf('Step: %d\n',it);
    u_new = u; % Copy previous solution

    % Compute diffusion using finite difference
    for i = 2:meshResX-1            % loop over the grid-x dir
        for j = 2:meshResY-1        % loop over the grid-y dir
            u_new(i, j) = u(i, j) + dt * nu * ...
                ((u(i+1, j) - 2*u(i, j) + u(i-1, j)) / meshsize^2 + ...
                (u(i, j+1) - 2*u(i, j) + u(i, j-1)) / meshsize^2);
        end
    end

    % Apply boundary conditions (Dirichlet - fixed zero at boundaries)
    u_new(1, :) = 0; u_new(end, :) = 0; % Left and right boundaries
    u_new(:, 1) = 0; u_new(:, end) = 0; % Top and bottom boundaries

    % Update solution
    u = u_new;
    save(sprintf('dat/step_%d.mat', it), 'X', 'Y', 'u', 'it', 't');
end

return
