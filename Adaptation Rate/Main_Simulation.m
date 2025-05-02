clear all
close all

addpath('AuxiliaryFunctions'); % contains all functions used by the main functions and scripts

%% Initialization ===========================================================================================
[initialPopulation, modelParameters, simulationParameters, discretizationParamaters, epsilon, eta] = initializeSimulation;

storageTimes = simulationParameters.times;
Dt = discretizationParamaters.Dt;
Dx = discretizationParamaters.Dx;
t = 0 : discretizationParamaters.Dt : simulationParameters.T;
x = simulationParameters.x_0 : Dx : simulationParameters.x_I;
I = length(x);

numVariables = length( fieldnames(initialPopulation) );   %number of variables 

population = struct('density', [], 'trait_mean', [], 'trait_variance', []);

population.density = zeros(I, length(storageTimes)); 
population.trait_mean = zeros(I, length(storageTimes)); 
population.trait_variance = zeros(I, length(storageTimes));

population.density(:,1) = initialPopulation.density; 
population.trait_mean(:,1) = initialPopulation.trait_mean; 
population.trait_variance(:,1) = initialPopulation.trait_variance;

%% Two-stage stabilizing correction ADI scheme [Hundsdorfer, App Num Math, 2000, p.222 =======================
systemSize = numVariables;
%---finding the location of the elements of solution matrices in U. 
%---varLocations(i,j) = linear index of the element at j_th row (j_th grid point) and i_th column (i_th variable) in solution matrix 
variableLocation = zeros(systemSize, I);
for i = 1:systemSize
    variableLocation(i,:) = [1 : systemSize : systemSize * I ] + (i-1);
end

solution = zeros(I,systemSize);
%---initializing solution equal to the initial values
solution(:, 1) = population.density(:,1);
solution(:, 2) = population.trait_mean(:,1);
solution(:, 3) = population.trait_variance(:,1);

odeSystemSize = numel(solution);
U = zeros(odeSystemSize,1);
for i = 1 : systemSize
    U(variableLocation(i,:)) = solution(:,i);
end

solution_star = zeros(I, systemSize);
storageCounter = 1;

figure, fig1 = axes;
% figure, fig2 = axes;
% figure, fig3 = axes;

randomAdaptationRate = zeros(I,1);
nonRandomAdaptationRate = zeros(I,1);

for k = 1 : length(t)-1
    
    plot(fig1, x, solution(:,1));
%     plot(fig2, x, solution(:,2));
%     plot(fig3, x, solution(:,3));
%  drawnow
    %tic
          
    %---Stage 1---------------------------------------------------------------------------------------------
    F1 = construct_F1(solution, variableLocation, Dx);
    [F2, A2] = construct_F2andA2(solution, variableLocation);
    F = F1 + F2;
    
    W0_star = U + Dt * F;
    
    A1 = construct_A1(solution, variableLocation, Dx);
    W1_star = ( speye(odeSystemSize) - eta * Dt * A1 ) \ ( W0_star - eta * Dt * A1 * U );
    
    U_star = ( speye(odeSystemSize) - eta * Dt * A2 ) \ ( W1_star - eta * Dt * A2 * U ); % U_star = W2_star
    
    for i = 1 : systemSize
        solution_star(:,i) =  U_star(variableLocation(i,:));
    end
      
    %---Stage 2---------------------------------------------------------------------------------------------
    F1 = construct_F1(solution_star, variableLocation, Dx);
    [F2, A2] = construct_F2andA2(solution_star, variableLocation);
    F_star = F1 + F2;
    
    W0 = U + 1/2 * Dt * ( F + F_star );
    
    A1 = construct_A1(solution_star, variableLocation, Dx);
    W1 = ( speye(odeSystemSize) - eta * Dt * A1 ) \ ( W0 - eta * Dt * A1 * U_star );
  
    U = ( speye(odeSystemSize) - eta * Dt * A2 ) \ ( W1 - eta * Dt * A2 * U_star ); % U = W2
    
    for i = 1 : systemSize
        solution(:,i) =  U(variableLocation(i,:));
    end
             
    %---solution storage-----------------------------------------------------
    if storageCounter + 1 <= length(storageTimes) 
        if abs( t(k+1) - storageTimes(storageCounter + 1) ) <= Dt/2
            storageCounter = storageCounter + 1;
            simulationParameters.times(storageCounter) = t(k+1); % actual storage times might be different that what set in simulation initialization (depending on the value of Dt)
            population.density(:,storageCounter) = solution(:,1);
            population.trait_mean(:,storageCounter) = solution(:,2);
            population.trait_variance(:,storageCounter) = solution(:,3);
        end 
    end
   % time = t(k+1) % Display current iteration time to monitor the progress
    %toc

    
    if sum(isnan(U)) % Means there is an error or the solution has diverged
        disp('Error! NaN entries are present in the solution')
        break
    end    
       
end 

%save('Data\sol_A_max_4_gradient_20.mat', 'population', 'modelParameters', 'simulationParameters', 'discretizationParamaters');


%save('Data\sol_A_max_4_V_20.mat', 'population', 'modelParameters', 'simulationParameters', 'discretizationParamaters');














