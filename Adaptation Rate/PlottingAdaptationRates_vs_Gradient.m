clear all
close all
addpath('AuxiliaryFunctions'); % contains all functions used by the main functions and scripts

numCurves = 1; %20;
%purple =  [0.4940, 0.1840, 0.5560];
purple =  [0.85, 0.6, 0.85];
%purple =  [0.6, 0.25, 0.7];

randomRate = zeros(numCurves,1);
nonRandomRate = zeros(numCurves,1);
totalRate = zeros(numCurves,1);
gradient = zeros(numCurves,1);

systemSize = 3;
for i = 1:numCurves
    path = strcat('Data\sol_A_max_10_gradient_', num2str(i) ,'.mat');
    load(path);
    
    N = population.density(:,end);
    Q = population.trait_mean(:,end);
    V = population.trait_variance(:,end);
    
    I = size(N,1);
    variableLocation = zeros(systemSize, I);
    for j = 1:systemSize
        variableLocation(j,:) = [1 : systemSize : systemSize * I ] + (j-1);
    end
    
    [random,  nonRandom] = calculateAdaptationRate([N Q V], variableLocation,  modelParameters, discretizationParamaters );
    random = random(variableLocation(2,:));
    nonRandom = nonRandom(variableLocation(2,:));
    total = random + nonRandom;
    
    peak = max(N);
    threshold = 0.49 * peak;  
    rangeLength = abs( floor(I/2) - dsearchn(N, threshold) );
    range = [floor(I/2): floor(I/2) + rangeLength];
    
    randomRate(i) = max(random(range));
    nonRandomRate(i) = max(nonRandom(range));
    totalRate(i) = max(total(range));
    
    gradient(i) = (modelParameters.Q_opt(end) - modelParameters.Q_opt(1)) / (simulationParameters.x_I-simulationParameters.x_0);
   
end


%%% Compute at least 4 curves and store them in the Data folder to be able
%%% to plot the graphs here. For the results shown in the paper, 20 curves
%%% were computed 

% fittedCurve = fit(gradient, totalRate,'exp2');
% 
% figure,
% hold on
% h = plot(fittedCurve, gradient, totalRate,'*');
% set(h, 'Color', purple);
% set(h, 'LineWidth', 1.25);
% hold off
% %grid(h, 'off')
% ylim([0, 0.25])
% xlabel('Gradient of the Optimal Trait  $[\mathtt{Q} / \mathtt{X}]$','Interpreter','latex','FontSize', 12);
% ylabel('Rate of Change in $q$ $[\mathtt{X}/\mathtt{T}]$','Interpreter','latex','FontSize', 12);
% grid on




