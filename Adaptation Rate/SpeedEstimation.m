clear all
close all

numCurves = 1; %20;
%purple =  [0.4940, 0.1840, 0.5560];
purple =  [0.85, 0.6, 0.85];

invasionSpeed = zeros(numCurves,1);
gradient = zeros(numCurves,1);
amplitude = zeros(numCurves,1);
maxVariance = zeros(numCurves,1);
for i = 1:numCurves
    path = strcat('Data\sol_A_max_10_gradient_', num2str(i) ,'.mat');
    load(path);
    n = population.density;
    amplitude(i) = max(n(:,end));
    gradient(i) = (modelParameters.Q_opt(end) - modelParameters.Q_opt(1)) / (simulationParameters.x_I-simulationParameters.x_0);
    times = simulationParameters.times;
    x = simulationParameters.x_0 : discretizationParamaters.Dx : simulationParameters.x_I;

    n = n(floor(size(n,1)/2):end, :); % only keep half of the symmetric density
    x = x(end-size(n,1)+1:end);
    
    threshold = 0.25 * amplitude(i);
    edgePoints = zeros(1,length(times));
    
    for j = 1:length(times)
        edge_index_x = dsearchn(n(:,j), threshold);
        edgePoints(j)= x(edge_index_x);
    end
    
    noThresholdCrossing = (edgePoints == simulationParameters.x_0);
    edgePoints(noThresholdCrossing) = []; % remove no threshold crossings
    times(noThresholdCrossing) = [];
    
    len = length(edgePoints);
    edgePoints(1:floor(len/5)) = []; %remove transient 20 percent at the begining
    times(1:floor(len/5)) = [];
    edgePoints(end-floor(len/10) : end) = []; % remove nonlinear 10 percent at the end
    times(end-floor(len/10) : end) = [];
    
    linearCoefficients = regress(edgePoints', [ones(length(times),1) times']); % fit a line
    invasionSpeed(i) = linearCoefficients(2);

    v = population.trait_variance;
    maxVariance(i) = v(ceil(length(v)/2), end); % maximum of v occurs at the center
end


%%% Compute at least 4 curves and store them in the Data folder to be able
%%% to plot the graphs here. For the results shown in the paper, 20 curves
%%% were computed 

figure,
plot(gradient, invasionSpeed, 'Color', purple, 'LineWidth', 1.5)
xlabel('Gradient of the Optimal Trait  $[\mathtt{Q} / \mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Invasion Wave Speed $[\mathtt{X}/\mathtt{T}]$','Interpreter','latex','FontSize', 12);
grid on
ylim([0,3.5]);

figure,
plot(gradient, amplitude, 'Color', purple, 'LineWidth', 1.5)
xlabel('Gradient of the Optimal Trait $[\mathtt{Q} / \mathtt{X}]$','Interpreter','latex','FontSize', 12);
ylabel('Invasion Wave Amplitude $[\mathtt{N}/\mathtt{X}]$','Interpreter','latex','FontSize', 12);
grid on
ylim([0,1.2]);



