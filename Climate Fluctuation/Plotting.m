clear all
close all

load('Results\sol_fluctuation_A_max_10_gradient_1_5')
fluctuationStartTime = 4; %[T] Q_opt starts fluctuating at this time


% Colors
darkBlue =  [0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];
yellow =  [0.9290, 0.6940, 0.1250];
purple =  [0.4940, 0.1840, 0.5560];
green =  [0.4660, 0.6740, 0.1880];
darkGreen =  [33/255, 186/255, 140/255];
lightBlue = [0.3010, 0.7450, 0.9330];
darkRed = [0.6350, 0.0780, 0.1840];
colors = [darkBlue; orange; yellow; purple; green; lightBlue; darkRed ];


colorMap1 = [linspace(160/256,256/256,64)' linspace(50/256,210/256,64)' linspace(0/256,145/256,64)']; % orange
%colorMap2 =[linspace(0.5,1,64)' linspace(0.2,0.9,64)' linspace(0.2,0.7,64)' ]; % pink
colorMap2 = [linspace(75/256,220/256,64)' linspace(110/256,256/256,64)' linspace(0/256,160/256,64)'];  % green
%colorMap2 =[linspace(0/256,195/256,64)' linspace(115/256,256/256,64)' linspace(100/256,220/256,64)' ]; % green


x_0 = simulationParameters.x_0;
x_I = simulationParameters.x_I;
Dx = discretizationParamaters.Dx;

x = x_0 : Dx : x_I;
I = length(x);
t = simulationParameters.times;

N = population.density;
Q = population.trait_mean;
V = population.trait_variance;

N(:, t<fluctuationStartTime) =[];
Q(:, t<fluctuationStartTime) =[];
V(:, t<fluctuationStartTime) =[];
t(:, t<fluctuationStartTime) =[];

numSamples = length(t);

figure, 
plot(x, N(:,[1:floor(numSamples/50):numSamples])','Color', orange, 'LineWidth', 0.5);
%plot(x, N(:,round(logspace(0,log10(numSamples), 20)) )','Color', orange, 'LineWidth', 0.5);
hold on
plot(x, N(:,1)','Color', orange, 'LineWidth', 1.5);
plot(x, N(:,201)','Color', darkRed, 'LineWidth', 1.5);
%plot(x, N(:,206)','Color', darkBlue, 'LineWidth', 1.5);
hold off
xlabel('Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Population Density $[\mathtt{N}/\mathtt{X}]$','Interpreter','latex','FontSize', 12);
%ylim = [0,1.2];
%colormap(fig_pop1_density, colorMap1);

figure, 
plot(x, Q(:,[1:floor(numSamples/20):numSamples])','Color', orange, 'LineWidth', 0.5);
hold on
plot(x, Q(:,1)','Color', orange, 'LineWidth', 1.5);
plot(x, Q(:,201)','Color', darkRed, 'LineWidth', 1.5);
plot(x, modelParameters.Q_opt','k', 'LineWidth', 1);
%plot(x, 2.5 + modelParameters.Q_opt','k--', 'LineWidth', 1);
hold off
xlabel('Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Trait Mean $[\mathtt{Q}]$','Interpreter','latex','FontSize', 12);
%ylim = [0,1.2];
%colormap(fig_pop1_density, colorMap1);

figure, 
plot(x, V(:,[1:floor(numSamples/20):numSamples])', 'Color', orange, 'LineWidth', 0.5);
hold on
plot(x, V(:,1)','Color', orange, 'LineWidth', 1.5);
plot(x, V(:,201)','Color', darkRed, 'LineWidth', 1.5);
hold off
xlabel('Space $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Trait Variance $[\mathtt{Q}^2]$','Interpreter','latex','FontSize', 12);
%ylim = [0,1.2];
%colormap(fig_pop1_density, colorMap1);



