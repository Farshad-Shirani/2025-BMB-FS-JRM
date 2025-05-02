function [initialPopulation, modelParameters, simulationParameters, discretizationParamaters, epsilon, eta] = initializeSimulation

epsilon = 0; %1e-14; % threshol for smallest value of population density (to avoid singularity due to 1/N term) 

%---simulation parameters---------------------------------------------------------------------------
T = 40; %final time
storageTimes = 0 : 0.1 : T; % solution storage times. Solutins are stored at closest time samples depending on Dt
x_0 = -50;    x_I = 50; % x1 range
simulationParameters = struct('T', T, 'times', storageTimes, 'x_0', x_0, 'x_I', x_I);

%---discretization parameters-----------------------------------------------------------------------
Dt = 0.002;  % time steps
Dx = (x_I - x_0)/ 1000; %x1 mesh   
discretizationParamaters = struct('Dt', Dt, 'Dx', Dx);

x = x_0 : Dx : x_I;
I = length(x);

%---model parameters--------------------------------------------------------------------------------
global D V_s V_u U kappa K Q_opt R A_max dQ_tilde
D = 1;  % diffusion matrix 
V_s = 1/0.2;    % 1/measure of the strength of stabilizing selection (V_s = 1/S, where S defined in Table 1)
V_u = 4;    % variance of the within phenotipic-resource utility curve (V_u = V where V is defined in Table 1)
U = 0.02;  % rate of increase in trait variance due to mutation
kappa = 0;  % population impact factor
K = 1 * ones(I, 1);  % constant carrying capacity
R = 2 * ones(I, 1);  % intrinsic rate of increase for optimally adapted individuals
dQ_opt = 0.2; % optimum trait gradient dQ
Q_0 = 10; % optimum trait at x_0
Q_opt = Q_0 + dQ_opt  * (x - x(1))';    % optimum trait value  

A_max = 10 * ones(I, 1); % maximum velocity adjustment for optimal advective dispersal (coefficient of the optimal dispersal term)

%---building perception of dQ_opt
dQ_max = 1; % maximum perceived ||dQ||
delta = 2; %width of the smooth cut-off function used at vicinity of the boundary
chi = @(x, a) ( exp( (x-a).^2 ./ ((x-a).^2 - delta^2 - eps) ) .* (1 - heaviside(abs(x-a) - delta) ) )'; % smooth cut-off function
dQ = gradient(Q_opt, Dx); 
dQ_norm = sqrt(sum(dQ.^2, 2));
dQ_cut = dQ - chi(x, x_0) * dQ(1) - chi(x, x_I) * dQ(end); % cutting the values of dQ at the vicinity of the boundaries
dQ_tilde = ( dQ_max ./ ( dQ_max + dQ_norm ) ) .* dQ_cut;

modelParameters = struct('D', D, 'V_s', V_s, 'V_u', V_u, 'U', U, 'kappa', kappa, 'K', K, 'Q_opt', Q_opt, 'R', R, 'A_max', A_max, 'dQ_tilde', dQ_tilde);

%---solver parameters-------------------------------------------------------------------------------
eta = 1/2;

%% Initial values ----------------------------------------------------------------------------------
initialPopulation = struct('density', [], 'trait_mean', [], 'trait_variance', []);
center = 0; % population center

%---For calculating nominal solution---------------------------------------------------------------- 
radius = sqrt(2); % population radius in each direction 
initialPopulation.density = 0.5 * sech( abs(x - center) / radius )' + epsilon;
center_index_x = dsearchn(x', center); 
initialPopulation.trait_mean = 0.6 * ( Q_opt - Q_opt(center_index_x) ) + Q_opt(center_index_x);
initialPopulation.trait_variance = 1 * ones(I, 1);


 




