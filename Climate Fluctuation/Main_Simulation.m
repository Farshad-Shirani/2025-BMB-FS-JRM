%---This is the extension of Case, 2000 model in 2D with variable variance.
clear all
close all

addpath('AuxiliaryFunctions'); % contains all functions used by the main functions and scripts


fluctuationDuration = [1, 1]; % [high pulse duration, low pulse duration]
%fluctuationDuration = [0.25, 0.25]; % [high pulse duration, low pulse duration]
fluctuationAmplitude = 5;
fluctuationStartTime = 4; %15 %[T] Q_opt starts fluctuating at this time

%% Initialization ===========================================================================================
[initialPopulation, modelParameters, simulationParameters, discretizationParamaters, epsilon, eta] = initializeSimulation;
global Q_opt

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
%---finding the location of the elements of solution matrices in U. Column no. in varLocations = linear index of the element in solution matrix 
variableLocation = zeros(systemSize, I);
for i = 1:systemSize
    variableLocation(i,:) = [1 : systemSize : systemSize * I] + (i-1);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilonFlag = zeros(length(t),1);

fluctuationSign = 1; % gives the sign of next alternating fluctuation in Q_opt
%Q_opt = Q_opt - fluctuationSign * fluctuationAmplitude / 2; % to adiust the first fluctuation amplitude
t0 = fluctuationStartTime;

for k = 1 : length(t)-1
    
    plot(fig1, x, solution(:,1));
%     plot(fig2, x, solution(:,2));
%     plot(fig3, x, solution(:,3));
%    drawnow
    %tic
    
    %---Fluctuating Q_opt-----------------------------------------
    if t(k) >= fluctuationStartTime
        if mod( (t(k)-t0), fluctuationDuration(1 + (1 + fluctuationSign)/2 ) ) < (0.98*Dt) % 0.98 is used because some times it happens that t(k)-t(0) is almost an espsion smaller than Dt in two consecutive steps
            Q_opt = Q_opt + fluctuationSign * fluctuationAmplitude;
            fluctuationSign = -fluctuationSign;
            t0 = t(k);
        end
    end
          
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
%    time = t(k+1) % Display current iteration time to monitor the progress
    %toc

    
    if sum(isnan(U)) % Means there is an error or the solution has diverged
        disp('Error! NaN entries are present in the solution')
        break
    end    
       
end 

%save('Results\sol_fluctuation_A_max_4_gradient_1_5.mat', 'population', 'modelParameters', 'simulationParameters', 'discretizationParamaters');
%save('Results\sol_fluctuation_A_max_4_gradient_1_5_Amp_9.mat', 'population', 'modelParameters', 'simulationParameters', 'discretizationParamaters');

%save('Results\fluctuation_vs_gradient\sol_fluctuation_A_max_10_gradient_1_8.mat', 'population', 'modelParameters', 'simulationParameters', 'discretizationParamaters');

%save('Results\fluctuation_vs_amplitude\sol_fluctuation_A_max_4_amplitude_4.mat', 'population', 'modelParameters', 'simulationParameters', 'discretizationParamaters');

%save('Results\fluctuation_vs_selection\sol_fluctuation_A_max_0_selection_0_2.mat', 'population', 'modelParameters', 'simulationParameters', 'discretizationParamaters');

%save('Results\fluctuation_vs_growthRate\sol_fluctuation_A_max_0_growthRate_2.mat', 'population', 'modelParameters', 'simulationParameters', 'discretizationParamaters');























