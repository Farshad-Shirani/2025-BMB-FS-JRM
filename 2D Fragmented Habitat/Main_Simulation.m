clear all
close all

addpath('AuxiliaryFunctions'); % contains all functions used by the main functions and scripts


%% Initialization ===========================================================================================
[initialPopulation, modelParameters, simulationParameters, discretizationParamaters, epsilon, eta] = initializeSimulation;

storageTimes = simulationParameters.times;
Dt = discretizationParamaters.Dt;
Dx1 = discretizationParamaters.Dx1;
Dx2 = discretizationParamaters.Dx2;
t = 0 : discretizationParamaters.Dt : simulationParameters.T;
x1 = simulationParameters.x1_0 : Dx1 : simulationParameters.x1_I;
x2 = simulationParameters.x2_0 : Dx2 : simulationParameters.x2_J;
I = length(x1);
J = length(x2);

if strcmp( simulationParameters.boundaryCondition_x2, 'periodic' )
    J_hat = J-1; % if periodic boundary condition is used in x2-direction
else
    J_hat = J; % if reflecting boundary condition is used in x2-direction
end

numVariables = length( fieldnames(initialPopulation) );   %number of variables 

population = struct('density', [], 'trait_mean', [], 'trait_variance', []);

population.density = zeros(I, J_hat, length(storageTimes)); 
population.trait_mean = zeros(I, J_hat, length(storageTimes)); 
population.trait_variance = zeros(I, J_hat, length(storageTimes));

population.density(:,:,1) = initialPopulation.density; 
population.trait_mean(:,:,1) = initialPopulation.trait_mean; 
population.trait_variance(:,:,1) = initialPopulation.trait_variance;

%% Two-stage stabilizing correction ADI scheme [Hundsdorfer, App Num Math, 2000, p.222 =======================
systemSize = numVariables;
%---finding the location of the elements of solution matrices in U. Column no. in varLocations = linear index of the element in solution matrix 
variableLocation = zeros(systemSize, I*J_hat);
for i = 1:systemSize
    variableLocation(i,:) = [1 : systemSize : systemSize * I * J_hat] + (i-1);
end

solution = zeros(I, J_hat, systemSize);
%---initializing solution equal to the initial values
solution(:, :, 1) = population.density(:,:,1);
solution(:, :, 2) = population.trait_mean(:,:,1);
solution(:, :, 3) = population.trait_variance(:,:,1);

solution_permuted = permute(solution, [2 1 3]); % transposes each solution matrix so that linear indexing orders the array in x2 direction

odeSystemSize = numel(solution);
U = zeros(odeSystemSize,1);
U_permuted = zeros(odeSystemSize,1); % is used at each stage when we alternate direction
for i = 1 : systemSize
    U(variableLocation(i,:)) = reshape(solution(:,:,i), [], 1);
    U_permuted(variableLocation(i,:)) = reshape(solution_permuted(:,:,i), [], 1);
end

solution_star = zeros(I, J_hat, systemSize);
U_star_permuted = zeros(odeSystemSize,1);
storageCounter = 1;

figure, fig1 = axes;
figure, fig2 = axes;
figure, fig3 = axes;
[X1, X2] = meshgrid(x1(1:I), x2(1:J_hat));

for k = 1 : length(t)-1
    
%     surf(fig1, X1, X2, solution(:,:,1)');
%     surf(fig2, X1, X2, solution(:,:,2)');
%     surf(fig3, X1, X2, solution(:,:,3)');
%     drawnow
    %tic
          
    %---Stage 1---------------------------------------------------------------------------------------------
    F1 = construct_F1(solution, variableLocation, Dx1);
    F2 = construct_F2(solution, variableLocation, Dx2);
    [F3, A3] = construct_F3andA3(solution, variableLocation);
    F = F1 + F2 + F3;
    
    W0_star = U + Dt * F;
    
    A1 = construct_A1(solution, variableLocation, Dx1);
    W1_star = ( speye(odeSystemSize) - eta * Dt * A1 ) \ ( W0_star - eta * Dt * A1 * U );
    
    %---alternating direction toward x2
    for i = 1 : systemSize
        buffer = W1_star(variableLocation(i,:));
        buffer = reshape(buffer, I, J_hat);
        buffer = buffer'; 
        W1_star(variableLocation(i,:)) = buffer(:);
    end
    
    A2 = construct_A2(solution_permuted, variableLocation, Dx2);
    W2_star = ( speye(odeSystemSize) - eta * Dt * A2 ) \ ( W1_star - eta * Dt * A2 * U_permuted );
    
    %---alternating direction back toward x1
    for i = 1 : systemSize
        buffer = W2_star(variableLocation(i,:));
        buffer = reshape(buffer, J_hat, I);
        buffer = buffer'; 
        W2_star(variableLocation(i,:)) = buffer(:);
    end
    
    U_star = ( speye(odeSystemSize) - eta * Dt * A3 ) \ ( W2_star - eta * Dt * A3 * U ); % U_star = W3_star
    
    for i = 1 : systemSize
        solution_star(:,:,i) =  reshape( U_star(variableLocation(i,:)), I, J_hat);
    end
    
    solution_star_permuted = permute(solution_star, [2 1 3]);
    for i = 1 : systemSize
        U_star_permuted(variableLocation(i,:)) = reshape(solution_star_permuted(:,:,i), [], 1);
    end
    
    %---Stage 2---------------------------------------------------------------------------------------------
    F1 = construct_F1(solution_star, variableLocation, Dx1);
    F2 = construct_F2(solution_star, variableLocation, Dx2);
    [F3, A3] = construct_F3andA3(solution_star, variableLocation);
    F_star = F1 + F2 + F3;
    
    W0 = U + 1/2 * Dt * ( F + F_star );
    
    A1 = construct_A1(solution_star, variableLocation, Dx1);
    W1 = ( speye(odeSystemSize) - eta * Dt * A1 ) \ ( W0 - eta * Dt * A1 * U_star );
    
    %---alternating direction toward x2
    for i = 1 : systemSize
        buffer = W1(variableLocation(i,:));
        buffer = reshape(buffer, I, J_hat);
        buffer = buffer'; 
        W1(variableLocation(i,:)) = buffer(:);
    end
    
    A2 = construct_A2(solution_star_permuted, variableLocation, Dx2);
    W2 = ( speye(odeSystemSize) - eta * Dt * A2 ) \ ( W1 - eta * Dt * A2 * U_star_permuted );
    
    %---alternating direction back toward x1
    for i = 1 : systemSize
        buffer = W2(variableLocation(i,:));
        buffer = reshape(buffer, J_hat, I);
        buffer = buffer'; 
        W2(variableLocation(i,:)) = buffer(:);
    end
    
    U = ( speye(odeSystemSize) - eta * Dt * A3 ) \ ( W2 - eta * Dt * A3 * U_star ); % U = W3
    
    for i = 1 : systemSize
        solution(:,:,i) =  reshape( U(variableLocation(i,:)), I, J_hat);
    end
        
    solution_permuted = permute(solution, [2 1 3]);
    for i = 1 : systemSize
        U_permuted(variableLocation(i,:)) = reshape(solution_permuted(:,:,i), [], 1);
    end
    
    %---solution storage-----------------------------------------------------
    if storageCounter + 1 <= length(storageTimes) 
        if abs( t(k+1) - storageTimes(storageCounter + 1) ) <= Dt/2
            storageCounter = storageCounter + 1;
            simulationParameters.times(storageCounter) = t(k+1); % actual storage times might be different that what set in simulation initialization (depending on the value of Dt)
            population.density(:,:,storageCounter) = solution(:,:,1);
            population.trait_mean(:,:,storageCounter) = solution(:,:,2);
            population.trait_variance(:,:,storageCounter) = solution(:,:,3);
        end 
    end
    time = t(k+1) % Display current iteration time to monitor the progress
    %toc

    
    if sum(isnan(U)) % Means there is an error or the solution has diverged
        disp('Error! NaN entries are present in the solution')
        break
    end    
       
end 
     

%save('Results\sol_patchyHabitat_gradient_1_A_0_D_4.mat', 'population', 'modelParameters', 'simulationParameters', 'discretizationParamaters');

%save('Results\sol_continuousHabitat_gradient_1_A_10.mat', 'population', 'modelParameters', 'simulationParameters', 'discretizationParamaters');




























