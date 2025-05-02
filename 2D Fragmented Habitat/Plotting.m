clear all
close all

addpath('AuxiliaryFunctions'); % contains all functions used by the main functions and scripts

load('Results\sol_patchyHabitat_gradient_1_A_10.mat')

colorMap1 = [linspace(160/255,255/255,64)' linspace(40/255,231/255,64)' linspace(5/256,209/255,64)' ]; % orange


x1_0 = simulationParameters.x1_0;
x1_I = simulationParameters.x1_I;
x2_0 = simulationParameters.x2_0;
x2_J = simulationParameters.x2_J;
Dx1 = discretizationParamaters.Dx1;
Dx2 = discretizationParamaters.Dx2;

x1 = x1_0 : Dx1 : x1_I;
x2 = x2_0 : Dx2 : x2_J;
I = length(x1);
J = length(x2);

[X1, X2] = meshgrid(x1(1:I), x2(1:J-1));
N = population.density;
Q = population.trait_mean;
V = population.trait_variance;
Q_optimum = modelParameters.Q_opt;

%% Plot sample frames ==================================================================
k = 21; % frame number

figure, fig_pop1_density = axes;
n1_plotHandle = surf(fig_pop1_density, X1, X2, N(:,:,k)');
xlabel(fig_pop1_density, '$x_1$ $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_pop1_density, '$x_2$ $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
zlabel(fig_pop1_density, '$Population Density$ $[\mathtt{N} / \mathtt{X}]$','Interpreter','latex','FontSize', 12);
set(n1_plotHandle,'EdgeColor', 'none'); % comment this line for 3D plot
% caxis(fig_pop1_density, [min(N(:)),max(N(:))] )
% colorbar
% axis off % comment this line for 3D plot
% axis image % comment this line for 3D plot
colormap(fig_pop1_density, colorMap1);


figure, fig_pop1_trait_mean = axes;
Q_difference = Q(:,:,k) - Q_optimum;
q1_plotHandle = surf(fig_pop1_trait_mean, X1, X2, Q_difference');
reduceEdgeLines(q1_plotHandle, 75, 75)
xlabel(fig_pop1_trait_mean, '$x_1$ $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_pop1_trait_mean, '$x_2$ $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
zlabel(fig_pop1_trait_mean, '$q - \mathrm{Q}$ $[\mathtt{Q}]$','Interpreter','latex','FontSize', 12);
set(q1_plotHandle,'EdgeColor', 'none');
%axis(fig_pop1_trait_mean, [min(Q_difference(:)), max(Q_difference(:))] )
%colorbar
% axis off
% axis image
colormap(fig_pop1_trait_mean, colorMap1);

figure, fig_pop1_trait_variance = axes;
v1_plotHandle = surf(fig_pop1_trait_variance, X1, X2, V(:,:,k)');
reduceEdgeLines(v1_plotHandle, 75, 75)
xlabel(fig_pop1_trait_variance, '$x_1$ $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel(fig_pop1_trait_variance, '$x_2$ $[\mathtt{X}]$','Interpreter','latex','FontSize', 12);
zlabel(fig_pop1_trait_variance, '$v$ $[\mathtt{Q}^2]$','Interpreter','latex','FontSize', 12);
set(v1_plotHandle,'EdgeColor', 'none');
%axis(fig_pop1_trait_variance, [min(V(:)), max(V(:))] )
%colorbar
% axis off
% axis image
colormap(fig_pop1_trait_variance, colorMap1);








