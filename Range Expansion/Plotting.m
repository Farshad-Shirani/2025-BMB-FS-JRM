clear all
close all

addpath('AuxiliaryFunctions'); % contains all functions used by the main functions and scripts

load('Results\sol_A_max_10.mat')

numSamples = length(simulationParameters.times);

curveIncrement = floor(numSamples/20); % increments between plotted curves
highlightedCurve = 81; %index of the highlighted curve

% Colors
darkBlue =  [0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];
yellow =  [0.9290, 0.6940, 0.1250];
purple =  [0.4940, 0.1840, 0.5560];
green =  [0.4660, 0.6740, 0.1880];
darkGreen =  [33/255, 186/255, 140/255];
lightBlue = [0.3010, 0.7450, 0.9330];
darkRed = [0.6350, 0.0780, 0.1840];
lightRed = [247/255, 210/255, 217/255];
lightOrange = [1, 0.85, 0.8];
colors = [darkBlue; orange; yellow; purple; green; lightBlue; darkRed ];


x_0 = simulationParameters.x_0;
x_I = simulationParameters.x_I;
Dx = discretizationParamaters.Dx;

x = x_0 : Dx : x_I;
I = length(x);

N = population.density;
Q = population.trait_mean;
V = population.trait_variance;

figure, 
plot(x, N(:,[1:curveIncrement:numSamples])','Color', orange, 'LineWidth', 0.5);
%plot(x, N(:,round(logspace(0,log10(numSamples), 20)) )','Color', orange, 'LineWidth', 0.5);
hold on
plot(x, N(:,1)','Color', orange, 'LineWidth', 1.75);
plot(x, N(:,highlightedCurve)','Color', darkRed, 'LineWidth', 1.75);
hold off
xlabel('Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Population Density $[\mathtt{N}/\mathtt{X}]$','Interpreter','latex','FontSize', 12);
%ylim = [0,1.2];
%colormap(fig_pop1_density, colorMap1);

figure, 
plot(x, Q(:,[1:curveIncrement:numSamples])','Color', orange, 'LineWidth', 0.5);
hold on
plot(x, Q(:,1)','Color', orange, 'LineWidth', 1.75);
plot(x, Q(:,highlightedCurve)','Color', darkRed, 'LineWidth', 1.75);
plot(x, modelParameters.Q_opt','k', 'LineWidth', 1.5);
%plot(x, 2.5 + modelParameters.Q_opt','k--', 'LineWidth', 1);
hold off
xlabel('Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Trait Mean $[\mathtt{Q}]$','Interpreter','latex','FontSize', 12);
%ylim = [0,1.2];
%colormap(fig_pop1_density, colorMap1);

figure, 
plot(x, V(:,[1:curveIncrement:numSamples])', 'Color', orange, 'LineWidth', 0.5);
hold on
plot(x, V(:,1)','Color', orange, 'LineWidth', 1.75);
plot(x, V(:,highlightedCurve)','Color', darkRed, 'LineWidth', 1.75);
hold off
xlabel('Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Trait Variance $[\mathtt{Q}^2]$','Interpreter','latex','FontSize', 12);
%ylim = [0,1.2];
%colormap(fig_pop1_density, colorMap1);


%% Calculate adptation rate componets ==============================================================
systemSize = 3;
I = size(N,1);
variableLocation = zeros(systemSize, I);
for i = 1:systemSize
    variableLocation(i,:) = [1 : systemSize : systemSize * I ] + (i-1);
end

figure, fig_random = axes;
figure, fig_nonRandom = axes;
figure, fig_total = axes;
hold(fig_random, 'on');
hold(fig_nonRandom, 'on');
hold(fig_total, 'on');
for i = 1 : 10 : size(N,2)
    [random,  nonRandom] = calculateAdaptationRate([N(:,i) Q(:,i) V(:,i)], variableLocation,  modelParameters, discretizationParamaters );
    randomAdaptationRate = random(variableLocation(2,:));
    nonRandomAdaptationRate = nonRandom(variableLocation(2,:));
    
    lineWidth = 0.5;
    lineColor_1 = lightOrange;
    lineColor_2 = orange;
    
    plot(fig_random, x, randomAdaptationRate','Color', lightOrange, 'LineWidth', lineWidth);
    plot(fig_nonRandom, x, nonRandomAdaptationRate','Color', lightOrange, 'LineWidth', lineWidth);
    plot(fig_total, x, randomAdaptationRate' + nonRandomAdaptationRate','Color', lightOrange, 'LineWidth', lineWidth);
    
    %---detecting edge-------------
    peak = max(N(:,i));
    threshold = 0.5 * peak;  
    rangeLength = abs( floor(I/2) - dsearchn(N(:,i), threshold) );
    range = [floor(I/2): floor(I/2) + rangeLength]; % right side of the range
    
    plot(fig_random, x(range), randomAdaptationRate(range)','Color', orange, 'LineWidth', lineWidth);
    plot(fig_nonRandom, x(range), nonRandomAdaptationRate(range)','Color', orange, 'LineWidth', lineWidth);
    plot(fig_total, x(range), randomAdaptationRate(range)' + nonRandomAdaptationRate(range)','Color', lineColor_2, 'LineWidth', lineWidth);
end
%---ploting the first and highlighted curves again so that they are shown in front
for i = 1 : 10 : size(N,2)
    [random,  nonRandom] = calculateAdaptationRate([N(:,i) Q(:,i) V(:,i)], variableLocation,  modelParameters, discretizationParamaters );
    randomAdaptationRate = random(variableLocation(2,:));
    nonRandomAdaptationRate = nonRandom(variableLocation(2,:));
    if i == 1
        lineWidth = 1.75;
        lineColor_1 = lightOrange;
        lineColor_2 = orange;
    elseif i == highlightedCurve 
        lineWidth = 1.75;
        lineColor_1 = lightRed;
        lineColor_2 = darkRed;
    else
        continue
    end
  
    plot(fig_random, x, randomAdaptationRate','Color', lineColor_1, 'LineWidth', lineWidth);
    plot(fig_nonRandom, x, nonRandomAdaptationRate','Color', lineColor_1, 'LineWidth', lineWidth);
    plot(fig_total, x, randomAdaptationRate' + nonRandomAdaptationRate','Color', lineColor_1, 'LineWidth', lineWidth);
    
    %---detecting edge-------------
    peak = max(N(:,i));
    threshold = 0.5 * peak;  
    rangeLength = abs( floor(I/2) - dsearchn(N(:,i), threshold) );
    range = [floor(I/2): floor(I/2) + rangeLength]; % right side of the range
    
    plot(fig_random, x(range), randomAdaptationRate(range)','Color', lineColor_2, 'LineWidth', lineWidth);
    plot(fig_nonRandom, x(range), nonRandomAdaptationRate(range)','Color', lineColor_2, 'LineWidth', lineWidth);
    plot(fig_total, x(range), randomAdaptationRate(range)' + nonRandomAdaptationRate(range)','Color', lineColor_2, 'LineWidth', lineWidth);
end

hold(fig_random, 'off');
hold(fig_nonRandom, 'off');
hold(fig_total, 'off');

xLimit = 40;
xlim(fig_random,[0, xLimit]);
xlim(fig_nonRandom,[0, xLimit]);
xlim(fig_total,[0, xLimit]);

xlabel(fig_random, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_random, 'Rate of Change in $q$ $[\mathtt{Q}/\mathtt{T}]$','Interpreter','latex','FontSize', 12);    
    
xlabel(fig_nonRandom, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_nonRandom, 'Rate of Change in $q$ $[\mathtt{Q}/\mathtt{T}]$','Interpreter','latex','FontSize', 12);    

xlabel(fig_total, 'Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_total, 'Rate of Change in $q$ $[\mathtt{Q}/\mathtt{T}]$','Interpreter','latex','FontSize', 12);    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    