%---calculating the fluctuations in populations density (max and min) at the steady state (after a
%long time). fluctuations are calculated at the center of the habitat
clear all
close all

%load('Results\fluctuation_vs_gradient\sol_fluctuation_A_max_10_gradient_1_8');
%load('Results\fluctuation_vs_amplitude\sol_fluctuation_A_max_4_amplitude_4');
%load('Results\fluctuation_vs_selection\sol_fluctuation_A_max_0_selection_0_2');
load('Results\fluctuation_vs_growthRate\sol_fluctuation_A_max_0_growthRate_2');


t = simulationParameters.times;
steadyStateInterval = [round(0.9*length(t)) : length(t)]; % last 10 percent of the simulation

N = population.density;
habitatLength = size(N,1);
N_center = N(round(habitatLength/2), :);

plot(t, N_center);

density_min = min( N_center(steadyStateInterval) )
density_max = max( N_center(steadyStateInterval) )