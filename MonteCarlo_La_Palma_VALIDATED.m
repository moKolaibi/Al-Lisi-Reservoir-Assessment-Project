% geothermal_montecarlo_la_palma_validated.m
% -----------------------------------------------------------
% Monte Carlo Simulation for Geothermal Power Plant Estimation
% Model Validation Against La Palma Study Parameters
%
% Author: Mugahed Kolaibi (moKolaibi)
% Date: 2025-09-26
%
% Description:
%   - Performs Monte Carlo sampling using La Palma (Canary Islands) field data
%   - Computes and visualizes total thermal energy, power harvested, and energy breakdown
%   - All results in MWe or MJ as described
%
% Usage:
%   1. Run this script to reproduce La Palma case validation.
%   2. See code comments for parameter distributions.
% -----------------------------------------------------------
% Number of simulations
num_simulations = 100000;

% Define fixed parameters
plant_life_years = 25; % Plant life (years)
plant_life_seconds = plant_life_years * 365 * 24 * 3600; % Plant life (seconds)
Cw = 4.18 * 1000; % Fluid specific heat (J/kg°C) - converted from kJ/kg°C

% Preallocate arrays for results
qR_values = zeros(num_simulations, 1);
qF_values = zeros(num_simulations, 1);
qT_values = zeros(num_simulations, 1);
Pe_values = zeros(num_simulations, 1);
Rf_values = zeros(num_simulations, 1);
eta_c_values = zeros(num_simulations, 1);
F_values = zeros(num_simulations, 1);

% Define distribution parameters
area_params = [6.5, 7.2, 13.7]; % Area (km^2)
thickness_params = [0.5, 0.9, 1.26]; % Thickness (km)
TR_uniform_params = [130, 200]; % Reservoir temperature (°C)
rho_r_params = [2000, 2300, 2600]; % Rock Density (kg/m^3)
phi_params = [0.01, 0.05, 0.1]; % Porosity (fraction)
Tref_params = [70, 85, 100]; % Reference temperature (°C)
Rf_params = [0, 10, 20]; % Recovery factor (%)
F_params = [80, 92.7, 95]; % Capacity factor (%)

% Monte Carlo Simulation
for i = 1:num_simulations
    % Sample from distributions
    area = randomTriangular(area_params) * 1e6; % Convert area from km^2 to m^2
    thickness = randomTriangular(thickness_params) * 1e3; % Convert thickness from km to m
    TR = randomUniform(TR_uniform_params);
    rho_r = randomTriangular(rho_r_params);
    phi = randomTriangular(phi_params);
    Tref = randomTriangular(Tref_params);
    Rf = randomTriangular(Rf_params) / 100; % Convert to fraction
    F = randomTriangular(F_params) / 100; % Convert to fraction
    
    % Calculate fluid density (kg/m^3)
    A = 1.16849 - 0.001477 * (TR + 273.15);
    rho_f = 1000 / (A + 3.05644e-6 * (TR + 273.15)^2);
    
    % Calculate rock specific heat (J/kg°C)
    Cr = (-4.418e-7 * TR^3 - 0.0008209 * TR^2 + 1.352 * TR + 994.2); % Convert from kJ/kg°C to J/kg°C
    
    % Calculate reservoir volume (m^3)
    V = area * thickness;
    
    % Calculate thermal energy stored in the rock and geofluid (J)
    qR = rho_r * Cr * V * (1 - phi) * (TR - Tref);
    qF = rho_f * Cw * V * phi * (TR - Tref);
    qT = qR + qF;
    
    % Convert qT from J to MJ
    qT = qT / 1e6;
    
    % Store values
    qR_values(i) = qR;
    qF_values(i) = qF;
    qT_values(i) = qT;
    Rf_values(i) = Rf;
    F_values(i) = F;
    
    % Calculate conversion efficiency (dependent on TR)
    eta_c = 0.0935 * TR - 2.3266; % Conversion efficiency in %
    eta_c = eta_c / 100; % Convert to fraction
    eta_c_values(i) = eta_c; % Store value
    
    % Calculate electrical power harvested (MWe)
    Pe = (qT * Rf * eta_c) / (F * plant_life_seconds);
    Pe_values(i) = Pe;
end

% Calculate statistics
mean_qR = mean(qR_values);
std_qR = std(qR_values);
min_qR = min(qR_values);
max_qR = max(qR_values);

mean_qF = mean(qF_values);
std_qF = std(qF_values);
min_qF = min(qF_values);
max_qF = max(qF_values);

mean_qT = mean(qT_values);
std_qT = std(qT_values);
min_qT = min(qT_values);
max_qT = max(qT_values);

mean_Pe = mean(Pe_values);
std_Pe = std(Pe_values);
min_Pe = min(Pe_values);
max_Pe = max(Pe_values);

mean_Rf = mean(Rf_values);
std_Rf = std(Rf_values);
min_Rf = min(Rf_values);
max_Rf = max(Rf_values);

mean_eta_c = mean(eta_c_values);
std_eta_c = std(eta_c_values);
min_eta_c = min(eta_c_values);
max_eta_c = max(eta_c_values);

mean_F = mean(F_values);
std_F = std(F_values);
min_F = min(F_values);
max_F = max(F_values);

% Display results
fprintf('Thermal Energy Stored in Rock (qR):\n');
fprintf('Mean: %.2e J\n', mean_qR);
fprintf('Standard Deviation: %.2e J\n', std_qR);
fprintf('Minimum: %.2e J\n', min_qR);
fprintf('Maximum: %.2e J\n\n', max_qR);

fprintf('Thermal Energy Stored in Geofluid (qF):\n');
fprintf('Mean: %.2e J\n', mean_qF);
fprintf('Standard Deviation: %.2e J\n', std_qF);
fprintf('Minimum: %.2e J\n', min_qF);
fprintf('Maximum: %.2e J\n\n', max_qF);

fprintf('Total Thermal Energy (qT):\n');
fprintf('Mean: %.2e MJ\n', mean_qT);
fprintf('Standard Deviation: %.2e MJ\n', std_qT);
fprintf('Minimum: %.2e MJ\n', min_qT);
fprintf('Maximum: %.2e MJ\n\n', max_qT);

fprintf('Electrical Power Harvested (Pe):\n');
fprintf('Mean: %.2f MWe\n', mean_Pe);
fprintf('Standard Deviation: %.2f MWe\n', std_Pe);
fprintf('Minimum: %.2f MWe\n', min_Pe);
fprintf('Maximum: %.2f MWe\n', max_Pe);

fprintf('Recovery Factor (Rf):\n');
fprintf('Mean: %.2f %%\n', mean_Rf * 100);
fprintf('Standard Deviation: %.2f %%\n', std_Rf * 100);
fprintf('Minimum: %.2f %%\n', min_Rf * 100);
fprintf('Maximum: %.2f %%\n\n', max_Rf * 100);

fprintf('Conversion Efficiency (eta_c):\n');
fprintf('Mean: %.2f %%\n', mean_eta_c * 100);
fprintf('Standard Deviation: %.2f %%\n', std_eta_c * 100);
fprintf('Minimum: %.2f %%\n', min_eta_c * 100);
fprintf('Maximum: %.2f %%\n\n', max_eta_c * 100);

fprintf('Plant Capacity Factor (F):\n');
fprintf('Mean: %.2f %%\n', mean_F * 100);
fprintf('Standard Deviation: %.2f %%\n', std_F * 100);
fprintf('Minimum: %.2f %%\n', min_F * 100);
fprintf('Maximum: %.2f %%\n\n', max_F * 100);

fprintf('Plant Lifetime:\n');
fprintf('%.2f seconds\n', plant_life_seconds);

% Save results for plotting in another script
save('simulation_results.mat', 'qR_values', 'qF_values', 'Pe_values', 'num_simulations');

% Helper function to generate random values from a triangular distribution
function val = randomTriangular(params)
    a = params(1);
    b = params(2);
    c = params(3);
    U = rand;
    F = (b - a) / (c - a);
    if U <= F
        val = a + sqrt(U * (b - a) * (c - a));
    else
        val = c - sqrt((1 - U) * (c - b) * (c - a));
    end
end

% Helper function to generate random values from a uniform distribution
function val = randomUniform(params)
    a = params(1);
    b = params(2);
    val = a + (b - a) * rand;
end

% Plotting the histogram with confidence intervals

% Define the number of bins for the histogram
num_bins = 50;

% Plot the frequency distribution of electrical power harvested
figure;
histogram(Pe_values, num_bins, 'Normalization', 'pdf');
xlabel('Electrical Power Harvested (GWe)');
ylabel('Probability Density');
title('Frequency Distribution of Electrical Power Harvested');
grid on;

% Set x-axis intervals
xticks = linspace(min(Pe_values), max(Pe_values), 11);
xticklabels = arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);

% Calculate 90% confidence interval
sorted_Pe_values = sort(Pe_values);
lower_bound = sorted_Pe_values(round(0.05 * num_simulations));
upper_bound = sorted_Pe_values(round(0.95 * num_simulations));

hold on;
xline(lower_bound, 'r--', 'Lower 90% CI');
xline(upper_bound, 'r--', 'Upper 90% CI');

% Display confidence interval on the plot
text(lower_bound, 0.8 * max(ylim), sprintf('%.2f', lower_bound), 'Color', 'red');
text(upper_bound, 0.8 * max(ylim), sprintf('%.2f', upper_bound), 'Color', 'red');

% Determine the most likely interval (mode of the distribution)
[counts, edges] = histcounts(Pe_values, num_bins, 'Normalization', 'pdf');
[~, max_idx] = max(counts);
most_likely_interval = [edges(max_idx), edges(max_idx + 1)];

% Highlight the most likely interval on the plot
patch([most_likely_interval(1), most_likely_interval(2), most_likely_interval(2), most_likely_interval(1)], ...
      [0, 0, max(counts), max(counts)], 'yellow', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

legend('Power Distribution', 'Lower 90% CI', 'Upper 90% CI', 'Most Likely Interval');
hold off;

% Plotting the cumulative frequency distribution

% Calculate the cumulative frequency
[cdf_values, cdf_edges] = histcounts(Pe_values, num_bins, 'Normalization', 'cdf');
cdf_centers = (cdf_edges(1:end-1) + cdf_edges(2:end)) / 2;

figure;
plot(cdf_centers, cdf_values, 'LineWidth', 2);
xlabel('Electrical Power Harvested (GWe)');
ylabel('Cumulative Frequency');
title('Cumulative Frequency Distribution of Electrical Power Harvested');
grid on;

% Calculate p10, p50, p90
p10 = prctile(Pe_values, 10);
p50 = prctile(Pe_values, 50);
p90 = prctile(Pe_values, 90);

hold on;
xline(p10, 'g--', 'P10');
xline(p50, 'b--', 'P50');
xline(p90, 'r--', 'P90');

% Display percentiles on the plot
text(p10, 0.1, sprintf('%.2f', p10), 'Color', 'green');
text(p50, 0.5, sprintf('%.2f', p50), 'Color', 'blue');
text(p90, 0.9, sprintf('%.2f', p90), 'Color', 'red');

legend('Cumulative Frequency', 'P10', 'P50', 'P90');
hold off;