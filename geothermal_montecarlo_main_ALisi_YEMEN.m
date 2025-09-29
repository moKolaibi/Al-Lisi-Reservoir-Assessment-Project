% geothermal_montecarlo_alisi_yemen.m
% -----------------------------------------------------------
% Monte Carlo Simulation for Geothermal Power Potential Estimation
% with Automated Convergence Testing
%
% Case Study: Al-Lisi Reservoir, Yemen
% Author: Mugahed Kolaibi (moKolaibi)
% Date: 2025-09-26
%
% Description:
%   - Runs a full Monte Carlo simulation for geothermal power capacity (MWe)
%   - Includes convergence diagnostics, statistics, and publication-ready plots
%   - Computes energy breakdown (rock/fluid/total), percentiles, and mode


clear; clc; rng(42);

%% Input Parameters
plant_life_years = 25;
plant_life_seconds = plant_life_years * 365 * 24 * 3600;
Cw    = 4180;    % Fluid specific heat [J/kg°C] (Fixed)
rho_f = 1000;    % Fluid density [kg/m³] (Fixed)
area_fixed = 200 * 1e6; % [m²] Fixed value

thick_min = 500; thick_mode = 650; thick_max = 800;              % [m] Triangular
Cr_min = 700; Cr_mode = 840; Cr_max = 1000;                      % [J/kg°C] Triangular
rho_r_min = 2500; rho_r_mode = 2700; rho_r_max = 3000;           % [kg/m³] Triangular
phi_min = 0.01; phi_max = 0.10;                                  % [fraction] Uniform
TR_min = 150; TR_mode = 200; TR_max = 250;                       % [°C] Triangular
Tref_min = 70; Tref_max = 100;                                   % [°C] Uniform
Rf_min = 0.03; Rf_mode = 0.11; Rf_max = 0.17;                    % [fraction] Triangular
F_min = 0.80; F_max = 0.95;                                      % [fraction] Uniform

num_bins = 40;

%% Convergence Testing Parameters
max_samples = 100000;  % Maximum number of samples
min_samples = 10000;    % Minimum samples before testing
check_interval = 5000;  % Check convergence every N samples
tolerance = 0.01;       % 1% tolerance for convergence
convergence_window = 3; % Number of consecutive checks needed

% Preallocate for maximum samples
PowerCapacity_MWe = zeros(max_samples,1);
qT_values = zeros(max_samples,1);
qF_values = zeros(max_samples,1);
qR_values = zeros(max_samples,1);

% Convergence tracking
mean_history = [];
converged = false;
convergence_count = 0;
final_sample_size = 0;

fprintf('Monte Carlo Geothermal Simulation - Started on 2025-09-26 12:40:35 UTC\n');
fprintf('User: moKolaibi\n');
fprintf('Starting Monte Carlo simulation with convergence testing...\n');
fprintf('Checking convergence every %d samples with %.2f%% tolerance\n', check_interval, tolerance*100);
fprintf('Maximum samples: %d, Minimum samples: %d\n', max_samples, min_samples);
fprintf('================================================\n');

%% Monte Carlo Loop with Convergence Testing
tic; % Start timing
for i = 1:max_samples
    thickness = randTriangular(thick_min, thick_mode, thick_max);         % [m]
    Cr        = randTriangular(Cr_min, Cr_mode, Cr_max);                 % [J/kg°C]
    rho_r     = randTriangular(rho_r_min, rho_r_mode, rho_r_max);        % [kg/m³]
    phi       = randInRange(phi_min, phi_max);                           % [fraction]
    TR        = randTriangular(TR_min, TR_mode, TR_max);                 % [°C]
    Tref      = randInRange(Tref_min, Tref_max);                         % [°C]
    Rf        = randTriangular(Rf_min, Rf_mode, Rf_max);                 % [fraction]
    F         = randInRange(F_min, F_max);                               % [fraction]

    V = area_fixed * thickness;
    qR = rho_r * Cr * V * (1 - phi) * (TR - Tref);      % [J]
    qF = rho_f * Cw * V * phi * (TR - Tref);            % [J]
    qT = qR + qF;                                       % [J]
    eta_c = max(0, 0.0935 * TR - 2.3266) / 100;         % Conversion efficiency (fraction)

    PowerCapacity_MWe(i) = (qT * Rf * eta_c) / (F * plant_life_seconds) / 1e6; % [MW]
    qT_values(i) = qT;
    qF_values(i) = qF;
    qR_values(i) = qR;

    % Check convergence
    if i >= min_samples && mod(i, check_interval) == 0
        current_mean = mean(PowerCapacity_MWe(1:i));
        mean_history = [mean_history, current_mean];
        
        if length(mean_history) >= 2
            relative_change = abs(mean_history(end) - mean_history(end-1)) / mean_history(end-1);
            
            if relative_change < tolerance
                convergence_count = convergence_count + 1;
                fprintf('Sample %d: Mean = %.2f MWe, Change = %.4f%% (Converged: %d/%d)\n', ...
                    i, current_mean, relative_change*100, convergence_count, convergence_window);
            else
                convergence_count = 0;
                fprintf('Sample %d: Mean = %.2f MWe, Change = %.4f%% (Not converged)\n', ...
                    i, current_mean, relative_change*100);
            end
            
            if convergence_count >= convergence_window
                fprintf('*** CONVERGENCE ACHIEVED at %d samples ***\n', i);
                converged = true;
                final_sample_size = i;
                break;
            end
        else
            fprintf('Sample %d: Mean = %.2f MWe (Initial check)\n', i, current_mean);
        end
    end
    
    % Progress indicator for large simulations
    if mod(i, 50000) == 0
        fprintf('Progress: %d samples completed...\n', i);
    end
end

simulation_time = toc;

if ~converged
    fprintf('*** Maximum samples (%d) reached without convergence ***\n', max_samples);
    final_sample_size = max_samples;
end

% Trim arrays to actual size used
PowerCapacity_MWe = PowerCapacity_MWe(1:final_sample_size);
qT_values = qT_values(1:final_sample_size);
qF_values = qF_values(1:final_sample_size);
qR_values = qR_values(1:final_sample_size);

num_samples = final_sample_size; % Update for rest of code

fprintf('\n================================================\n');
fprintf('SIMULATION COMPLETED\n');
fprintf('Final Results based on %d samples\n', num_samples);
if converged
    fprintf('Convergence Status: CONVERGED\n');
else
    fprintf('Convergence Status: NOT CONVERGED\n');
end
fprintf('Total Simulation Time: %.2f seconds\n', simulation_time);
fprintf('================================================\n');

%% Main Statistics (in MWe)
power_min = min(PowerCapacity_MWe);
power_max = max(PowerCapacity_MWe);
power_mode = modeFromHist(PowerCapacity_MWe, num_bins);
power_mean = mean(PowerCapacity_MWe);
power_std = std(PowerCapacity_MWe);
power_lb = prctile(PowerCapacity_MWe, 10);   % 10% CI
power_ub = prctile(PowerCapacity_MWe, 90);   % 90% CI
p10 = power_lb; % Use identical value for P10 everywhere
p50 = prctile(PowerCapacity_MWe, 50);
p90 = power_ub; % Use identical value for P90 everywhere

qT_min = min(qT_values); qT_max = max(qT_values);
qT_mode = modeFromHist(qT_values, num_bins);
qT_mean = mean(qT_values); qT_std = std(qT_values);
qT_lb = prctile(qT_values, 10);
qT_ub = prctile(qT_values, 90);
qT_p50 = prctile(qT_values, 50); % P50 for total energy

%% Energy Breakdown Statistics (Updated with P50)
qF_mean = mean(qF_values); qR_mean = mean(qR_values);
qF_mode = modeFromHist(qF_values, num_bins); qR_mode = modeFromHist(qR_values, num_bins);
qF_std = std(qF_values); qR_std = std(qR_values);

qF_p10 = prctile(qF_values, 10); qF_p50 = prctile(qF_values, 50); qF_p90 = prctile(qF_values, 90);
qR_p10 = prctile(qR_values, 10); qR_p50 = prctile(qR_values, 50); qR_p90 = prctile(qR_values, 90);

fluid_fraction_mean = qF_mean / qT_mean;
rock_fraction_mean = qR_mean / qT_mean;

%% Results Table
ResultsTable = table( ...
    [power_min; power_mean; power_mode; power_max; power_lb; power_ub; power_std], ...
    [qT_min; qT_mean; qT_mode; qT_max; qT_lb; qT_ub; qT_std], ...
    [p10; p50; p90; NaN; NaN; NaN; NaN], ...
    'VariableNames', {'PowerCapacity_MWe', 'TotalThermalEnergy_J', 'PowerPercentiles_MWe'}, ...
    'RowNames', {'Min', 'Mean', 'Mode', 'Max', '10% CI', '90% CI', 'Std/P10/P50/P90'});
disp('Main Results Table (in MWe):');
disp(ResultsTable);

%% Print Recommended Results (Geothermics style) - Updated with P50
fprintf('\n===== Geothermal Resource Breakdown =====\n');
fprintf('Mean energy stored in geofluid (qF): %.2e J\n', qF_mean);
fprintf('Mean energy stored in rock (qR):     %.2e J\n', qR_mean);
fprintf('Mean total energy (qT):              %.2e J\n', qT_mean);
fprintf('Fraction stored in fluid: %.2f%%\n', 100*fluid_fraction_mean);
fprintf('Fraction stored in rock:  %.2f%%\n', 100*rock_fraction_mean);
fprintf('Mode energy in fluid: %.2e J\n', qF_mode);
fprintf('Mode energy in rock:  %.2e J\n', qR_mode);
fprintf('Mode total energy:    %.2e J\n', qT_mode);
fprintf('P10 energy in fluid:  %.2e J\n', qF_p10);
fprintf('P50 energy in fluid:  %.2e J\n', qF_p50);
fprintf('P90 energy in fluid:  %.2e J\n', qF_p90);
fprintf('P10 energy in rock:   %.2e J\n', qR_p10);
fprintf('P50 energy in rock:   %.2e J\n', qR_p50);
fprintf('P90 energy in rock:   %.2e J\n', qR_p90);
fprintf('P10 total energy:     %.2e J\n', qT_lb);
fprintf('P50 total energy:     %.2e J\n', qT_p50);
fprintf('P90 total energy:     %.2e J\n', qT_ub);
if qF_mean > qR_mean
    fprintf('Geofluid stores more energy than rock on average.\n');
else
    fprintf('Rock stores more energy than geofluid on average.\n');
end
fprintf('=========================================\n');

%% Complete Energy Breakdown Table
fprintf('\n===== Complete Energy Breakdown Table =====\n');
fprintf('Statistic\tEnergy in Geofluid (qF)\tEnergy in Rock (qR)\tTotal Energy (qT)\tFraction Fluid (%%)\tFraction Rock (%%)\n');
fprintf('Mean\t\t%.2e J\t\t\t%.2e J\t\t\t%.2e J\t\t%.2f\t\t\t%.2f\n', qF_mean, qR_mean, qT_mean, 100*fluid_fraction_mean, 100*rock_fraction_mean);
fprintf('Mode\t\t%.2e J\t\t\t%.2e J\t\t\t%.2e J\t\t--\t\t\t--\n', qF_mode, qR_mode, qT_mode);
fprintf('P10\t\t\t%.2e J\t\t\t%.2e J\t\t\t%.2e J\t\t--\t\t\t--\n', qF_p10, qR_p10, qT_lb);
fprintf('P50\t\t\t%.2e J\t\t\t%.2e J\t\t\t%.2e J\t\t--\t\t\t--\n', qF_p50, qR_p50, qT_p50);
fprintf('P90\t\t\t%.2e J\t\t\t%.2e J\t\t\t%.2e J\t\t--\t\t\t--\n', qF_p90, qR_p90, qT_ub);
fprintf('============================================\n');

%% Convergence Summary
fprintf('\n===== Convergence Summary =====\n');
if converged
    fprintf('Convergence Status: CONVERGED\n');
else
    fprintf('Convergence Status: NOT CONVERGED\n');
end
fprintf('Final Sample Size: %d\n', num_samples);
fprintf('Convergence Tolerance: %.2f%%\n', tolerance*100);
fprintf('Consecutive Convergence Checks Required: %d\n', convergence_window);
fprintf('Total Simulation Time: %.2f seconds\n', simulation_time);
fprintf('Samples per second: %.0f\n', num_samples/simulation_time);
fprintf('===============================\n');

%% Plotting (PDF, CDF, Thermal) - all blue, all in MWe
geothermal_montecarlo_allblueplots_mwe('pdf', PowerCapacity_MWe, num_bins, [power_min, power_max, power_mode, power_std], [power_lb, power_ub], num_samples);
geothermal_montecarlo_allblueplots_mwe('cdf', PowerCapacity_MWe, num_bins, [p10, p50, p90], [], num_samples);
geothermal_montecarlo_allblueplots_mwe('thermal', qT_values, num_bins, [qT_min, qT_max, qT_mode, qT_mean, qT_std], [qT_lb, qT_ub], num_samples);

%% Plot Convergence History
if length(mean_history) > 1
    figure('Color','w','Name','Convergence History','NumberTitle','off');
    sample_points = min_samples:check_interval:(min_samples + (length(mean_history)-1)*check_interval);
    plot(sample_points, mean_history, 'b-', 'LineWidth', 2);
    xlabel('Number of Samples', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Mean Power Capacity (MWe)', 'FontSize', 14, 'FontWeight', 'bold');
    title('Monte Carlo Convergence History', 'FontSize', 16, 'FontWeight', 'bold');
    grid on;
    if converged
        hold on;
        plot(final_sample_size, mean_history(end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        text(final_sample_size, mean_history(end), sprintf('  Converged at %d samples', final_sample_size), ...
            'FontSize', 12, 'VerticalAlignment', 'middle');
    end
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');
    set(gcf, 'Position', [100 100 800 500]);
end

%% Helper Functions
function val = randInRange(a, b)
    val = a + (b-a)*rand;
end

function val = randTriangular(a, c, b)
    F = (c - a) / (b - a);
    u = rand;
    if u < F
        val = a + sqrt(u * (b - a) * (c - a));
    else
        val = b - sqrt((1 - u) * (b - a) * (b - c));
    end
end

function mode_val = modeFromHist(data, num_bins)
    [counts, edges] = histcounts(data, num_bins);
    [~, idx] = max(counts);
    mode_val = mean(edges(idx:idx+1));
end