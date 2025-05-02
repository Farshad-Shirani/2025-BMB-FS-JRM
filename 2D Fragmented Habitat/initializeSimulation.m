function [initialPopulation, modelParameters, simulationParameters, discretizationParamaters, epsilon, eta] = initializeSimulation

epsilon = 0;%1e-16; % threshol for smallest value of population density (to avoid singularity due to 1/N term) 

%---simulation parameters----------------------------------------------
T = 40; %final time
storageTimes = 0 : 2 : T; % solution storage times. Solutins are sored at closest time samples depending on Dt
x1_0 = -25;    x1_I = 25; % x1 range
x2_0 = -25;    x2_J = 25; % x2 range

global boundaryCondition_x2
boundaryCondition_x2 = 'periodic'; %boundary condition in x_2 direction
%boundaryCondition_x2 = 'reflecting'; %boundary condition in x_2 direction
simulationParameters = struct('T', T, 'times', storageTimes, 'x1_0', x1_0, 'x1_I', x1_I, 'x2_0', x2_0, 'x2_J', x2_J, 'boundaryCondition_x2', boundaryCondition_x2);

%---discretization parameters------------------------------------------
Dt = 0.01;  % time steps
Dx1 = (x1_I - x1_0)/ 400; %x1 mesh   
Dx2 = (x2_J - x2_0)/ 400; %x2 mesh  
discretizationParamaters = struct('Dt', Dt, 'Dx1', Dx1, 'Dx2', Dx2);

x1 = x1_0 : Dx1 : x1_I;
x2 = x2_0 : Dx2 : x2_J;
I = length(x1);
J = length(x2);

%---model parameters----------------------------------------------------
global D V_s V_u U kappa K Q_opt R A_max dQ_tilde
D = 1 * diag([1 1]);  % diffusion matrix 
V_s = 1/0.2;    % 1/measure of the strength of stabilizing selection (V_s = 1/S, where S defined in Table 1)
V_u = 4;    % variance of the within phenotipic-resource utility curve (V_u = V where V is defined in Table 1)
U = 0.02;  % rate of increase in trait variance due to mutation
kappa = 0;  % population impact factor

if strcmp( boundaryCondition_x2, 'periodic' )
    J_hat = J-1; % if periodic boundary condition is used in x2-direction
else
    J_hat = J; % if reflecting boundary condition is used in x2-direction
end


%K = 1 * ones(I, J_hat);  % carrying capacity
K = createPatchyHabitat(simulationParameters, discretizationParamaters)';

R = 2 * ones(I, J_hat);  % intrinsic rate of increase for optimally adapted individuals

dQ_opt = 1; % optimum trait gradient
Q_0 = 10; % optimum trait at x1_0
Q_opt = Q_0 + dQ_opt  * (x1 - x1(1))' * ones(1, J_hat);    % optimum trait value  
% [X1, X2] = meshgrid(x1(1:I), x2(1:J_hat));
% surf(X1,X2,Q_opt');

A_max = 0 * ones(I, J_hat); % maximum velocity adjustment for optimal advective dispersal (coefficient of the optimal dispersal term)

%---building perception of dQ_opt
dQ_max = 1; % maximum perceived ||dQ||
delta = 2; %width of the smooth cut-off function used at vicinity of the boundary
chi = @(x, a) ( exp( (x-a).^2 ./ ((x-a).^2 - delta^2 - eps) ) .* (1 - heaviside(abs(x-a) - delta) ) )'; % smooth cut-off function
[dQ_2, dQ_1] = gradient(Q_opt, Dx2, Dx1); 
dQ_norm = sqrt(dQ_1.^2 + dQ_2.^2);
%-cut values of dQ near the boundary. No need to cut if periodic boundary condition is used.
dQ_1_cut = dQ_1 - chi(x1, x1_0) * dQ_1(1,:)  - chi(x1, x1_I) * dQ_1(end,:); % cutting the values of dQ_1 at the vicinity of the boundaries
if strcmp( boundaryCondition_x2, 'periodic' )
    dQ_2_cut = dQ_2; % if periodic boundaru
else
    dQ_2_cut = dQ_2 - dQ_2(:,1) * chi(x2, x2_0)' - dQ_2(:,end) * chi(x2, x2_J)'; % cutting the values of dQ_2 at the vicinity of the boundaries
end
dQ_tilde(:,:,1) = ( dQ_max ./ ( dQ_max + dQ_norm ) ) .* dQ_1_cut; % dQ_tilde in x_1 direction
dQ_tilde(:,:,2) = ( dQ_max ./ ( dQ_max + dQ_norm ) ) .* dQ_2_cut; % dQ_tilde in x_2 direction


modelParameters = struct('D', D, 'V_s', V_s, 'V_u', V_u, 'U', U, 'kappa', kappa, 'K', K, 'Q_opt', Q_opt, 'R', R, 'A_max', A_max, 'dQ_tilde', dQ_tilde);

%---solver parameters----------------------------------------------------
eta = 1/2;

%---initial values--------------------------------------------------------
initialPopulation = struct('density', [], 'trait_mean', [], 'trait_variance', []);
[X1, X2] = meshgrid(x1(1:I), x2(1:J_hat));

center = [0 0]; % population center
radius = 2 * [1; 1]; % population radius in each direction 
initialPopulation.density = 1 * sech( sqrt( ( X1 - center(1) ).^2 / radius(1) + ( X2 - center(2) ).^2 / radius(2) ) )' + epsilon;

center_index_x1 = dsearchn(x1', center(1)); 
center_index_x2 = dsearchn(x2', center(2));
initialPopulation.trait_mean = 0.7 * ( Q_opt - Q_opt(center_index_x1, center_index_x2) ) + Q_opt(center_index_x1, center_index_x2);

initialPopulation.trait_variance = 1 * ones(I, J_hat);

end

