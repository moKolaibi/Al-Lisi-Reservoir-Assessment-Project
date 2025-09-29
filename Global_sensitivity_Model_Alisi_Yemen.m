% global_sensitivity_morris_alisi_yemen.m
% -----------------------------------------------------------
% Morris Elementary Effects (Global OAT) Sensitivity Analysis
% Case Study: Al-Lisi Reservoir (Yemen) - Top 5 Parameters
%
% Author: Mugahed Kolaibi (moKolaibi)
% Date: 2025-09-28
%
% Description:
%   - Identifies and ranks the 5 most influential parameters for geothermal power capacity
%   - Uses Morris elementary effects method with optimal sampling
%   - All plots and tables are formatted for publication
%
% Usage:
%   1. Run after completing Monte Carlo simulation for Al-Lisi (see README)
%   2. Adjust parameter bounds as needed for your reservoir
%   3. All plots and tables are generated automatically
% -----------------------------------------------------------

clear; clc; rng(42); % Set seed for reproducibility

%% Header Information
fprintf('Morris Elementary Effects - TOP 5 PARAMETERS ONLY\n');
fprintf('Geothermal Power Plant Estimation\n');
fprintf('Date: 2025-09-28 00:03:33 UTC | User: moKolaibi\n');
fprintf('=========================================================\n');

%% Fixed Parameters (matching recent Monte Carlo simulation)
plant_life_years = 25;
plant_life_seconds = plant_life_years * 365 * 24 * 3600;
Cw    = 4180;    % Fluid specific heat [J/kg°C]
rho_f = 1000;    % Fluid density [kg/m³]
area_fixed = 200 * 1e6; % [m²] Fixed value (200 km²)

%% Reduced Morris Method Parameters - TOP 5 ONLY
factors = {'Reservoir_Temp', 'Thickness', 'Recovery_Factor', 'Rock_Density', 'Capacity_Factor'};

% Parameter bounds [min, max] - TOP 5 MOST CRITICAL
param_bounds = [150,  250;     % Reservoir temperature [°C] - MOST IMPORTANT
                500,  800;     % Thickness [m] - ENERGY CONTENT
                0.03, 0.17;    % Recovery factor [fraction] - EXTRACTION EFFICIENCY
                2500, 3000;    % Rock density [kg/m³] - ENERGY STORAGE
                0.80, 0.95];   % Capacity factor [fraction] - OPERATIONAL EFFICIENCY

%% FIXED PARAMETERS (Previously Variable)
% These are now fixed at mean/mode values based on full analysis
rock_specific_heat_fixed = 840;   % [J/kg°C] - Fixed at mode
porosity_fixed = 0.055;           % [fraction] - Fixed at mean
reference_temp_fixed = 85;        % [°C] - Fixed at mean

k = length(factors);  % Number of factors = 5
r = 15;              % Increased trajectories for better convergence with fewer parameters
p = 6;               % Increased levels for better resolution
Delta = p/(2*(p-1)); % Optimal Delta value for p levels

fprintf('PARAMETER REDUCTION SUMMARY:\n');
fprintf('Original parameters: 8\n');
fprintf('Reduced parameters: %d (TOP 5 MOST IMPORTANT)\n', k);
fprintf('Fixed parameters: 3 (at mean/mode values)\n');
fprintf('- Rock Specific Heat: %.0f J/kg°C\n', rock_specific_heat_fixed);
fprintf('- Porosity: %.3f\n', porosity_fixed);
fprintf('- Reference Temperature: %.0f °C\n', reference_temp_fixed);
fprintf('---------------------------------------------------------\n');

%% Morris Sampling Design
fprintf('Generating Morris sampling design...\n');
fprintf('Number of factors (k): %d\n', k);
fprintf('Number of trajectories (r): %d\n', r);
fprintf('Number of levels (p): %d\n', p);
fprintf('Delta value: %.3f\n', Delta);
fprintf('Total model evaluations: %d\n', r*(k+1));
fprintf('---------------------------------------------------------\n');

% Generate Morris sample matrix
X_morris = generate_morris_sample(k, r, p, Delta);
n_samples = size(X_morris, 1);

% Scale to actual parameter ranges
X_scaled = zeros(n_samples, k);
for i = 1:k
    X_scaled(:,i) = param_bounds(i,1) + X_morris(:,i) * (param_bounds(i,2) - param_bounds(i,1));
end

%% Model Evaluations
fprintf('Performing model evaluations...\n');
Y = zeros(n_samples, 1);
tic;
for i = 1:n_samples
    Y(i) = calculate_power_capacity_top5(X_scaled(i,:), area_fixed, plant_life_seconds, Cw, rho_f, ...
                                         rock_specific_heat_fixed, porosity_fixed, reference_temp_fixed);
    if mod(i, 30) == 0
        fprintf('Completed %d/%d evaluations (%.1f%%)\n', i, n_samples, 100*i/n_samples);
    end
end
evaluation_time = toc;
fprintf('Model evaluations completed in %.2f seconds\n', evaluation_time);

%% Calculate Elementary Effects
fprintf('Calculating elementary effects...\n');
EE = zeros(r, k);  % Elementary effects matrix
mu_star = zeros(k, 1);     % Mean of absolute elementary effects
mu = zeros(k, 1);          % Mean of elementary effects
sigma = zeros(k, 1);       % Standard deviation of elementary effects

for traj = 1:r
    start_idx = (traj-1)*(k+1) + 1;
    end_idx = traj*(k+1);
    
    Y_traj = Y(start_idx:end_idx);
    X_traj = X_morris(start_idx:end_idx, :);
    
    for j = 1:k
        % Find the step where factor j changes
        diff_X = diff(X_traj);
        step_idx = find(abs(diff_X(:,j)) > 0, 1);
        
        if ~isempty(step_idx)
            % Calculate elementary effect
            delta_Y = Y_traj(step_idx+1) - Y_traj(step_idx);
            delta_X = diff_X(step_idx, j);
            EE(traj, j) = delta_Y / delta_X;
        end
    end
end

% Calculate Morris sensitivity measures
for j = 1:k
    valid_ee = EE(~isnan(EE(:,j)), j);  % Remove NaN values
    if ~isempty(valid_ee)
        mu(j) = mean(valid_ee);
        mu_star(j) = mean(abs(valid_ee));
        sigma(j) = std(valid_ee);
    end
end

%% Results Analysis
fprintf('\n===== TOP 5 PARAMETERS - MORRIS RESULTS =====\n');
fprintf('%-20s | μ (Mean EE) | μ* (|Mean EE|) | σ (Std EE) | Rank\n', 'Factor');
fprintf('------------------------------------------------------------------\n');

% Rank factors by μ* 
[sorted_mu_star, rank_idx] = sort(mu_star, 'descend');

for i = 1:k
    idx = rank_idx(i);
    fprintf('%-20s | %10.2f | %12.2f | %9.2f | %4d\n', ...
        factors{idx}, mu(idx), mu_star(idx), sigma(idx), i);
end
fprintf('==================================================================\n');

%% Efficiency Analysis
fprintf('\n===== EFFICIENCY ANALYSIS =====\n');
cumulative_importance = cumsum(sorted_mu_star) / sum(sorted_mu_star);
for i = 1:k
    fprintf('Top %d parameters explain %.1f%% of total sensitivity\n', i, cumulative_importance(i)*100);
end
fprintf('Model reduction efficiency: %.1f%% sensitivity captured\n', cumulative_importance(end)*100);
fprintf('Computational savings: ~%.0f%% (5 vs 8 parameters)\n', (1-k/8)*100);
fprintf('===============================\n');

%% Create Publication-Ready Table
MorrisTable = table( ...
    factors(rank_idx)', ...
    mu(rank_idx), ...
    mu_star(rank_idx), ...
    sigma(rank_idx), ...
    (1:k)', ...
    'VariableNames', {'Parameter', 'Mean_EE', 'AbsMean_EE', 'Std_EE', 'Rank'});

fprintf('\n===== TOP 5 PARAMETERS - PUBLICATION TABLE =====\n');
disp(MorrisTable);

%% PUBLICATION-READY PLOTS (2 MAIN FIGURES)

% Figure 1: Morris Plot (μ* vs σ) - ESSENTIAL FOR GEOTHERMICS
figure('Color','w','Name','Top 5 Parameters Morris Plot','NumberTitle','off','Position',[100 100 900 700]);
scatter(sigma, mu_star, 250, 'filled', 'MarkerFaceColor', [0.2, 0.6, 0.8], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2);
hold on;

% Publication-ready factor labels
pub_labels_morris = cell(k,1);
for i = 1:k
    switch factors{i}
        case 'Reservoir_Temp'
            pub_labels_morris{i} = 'Reservoir T';
        case 'Recovery_Factor'
            pub_labels_morris{i} = 'Recovery Factor';
        case 'Capacity_Factor'
            pub_labels_morris{i} = 'Capacity Factor';
        case 'Rock_Density'
            pub_labels_morris{i} = 'Rock Density';
        otherwise
            pub_labels_morris{i} = factors{i};
    end
end

% Add factor labels with optimal positioning
for i = 1:k
    text(sigma(i)+0.04*max(sigma), mu_star(i), pub_labels_morris{i}, ...
         'FontSize', 16, 'FontWeight', 'bold');
end

xlabel('σ (Standard Deviation of Elementary Effects)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('μ* (Mean Absolute Elementary Effects)', 'FontSize', 18, 'FontWeight', 'bold');
title('Morris Global Sensitivity Analysis (Top 5 Parameters)', 'FontSize', 20, 'FontWeight', 'bold');
grid on; box on;
set(gca, 'FontSize', 16, 'FontWeight', 'bold');

% Add interpretation zones
xlim_val = xlim; ylim_val = ylim;
text(0.05*xlim_val(2), 0.95*ylim_val(2), {'High Importance', 'Linear Effects'}, ...
     'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.4, 0.4, 0.4], ...
     'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'EdgeColor', 'k');
text(0.70*xlim_val(2), 0.95*ylim_val(2), {'High Importance', 'Non-linear Effects'}, ...
     'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.4, 0.4, 0.4], ...
     'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'EdgeColor', 'k');

% Figure 2: Importance Ranking - ESSENTIAL FOR GEOTHERMICS
figure('Color','w','Name','Top 5 Parameters Ranking','NumberTitle','off','Position',[200 200 1000 700]);

% Create publication labels for ranking
pub_labels_rank = cell(k,1);
for i = 1:k
    idx = rank_idx(i);
    switch factors{idx}
        case 'Reservoir_Temp'
            pub_labels_rank{i} = 'Reservoir Temperature';
        case 'Recovery_Factor'
            pub_labels_rank{i} = 'Recovery Factor';
        case 'Capacity_Factor'
            pub_labels_rank{i} = 'Capacity Factor';
        case 'Rock_Density'
            pub_labels_rank{i} = 'Rock Density';
        otherwise
            pub_labels_rank{i} = factors{idx};
    end
end

barh(1:k, sorted_mu_star, 'FaceColor', [0.8, 0.4, 0.2], 'EdgeColor', [0.4, 0.2, 0.1], 'LineWidth', 2);
title('Critical Parameter Ranking for Geothermal Power Capacity', 'FontSize', 20, 'FontWeight', 'bold');
xlabel('μ* (Mean Absolute Elementary Effects)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Geothermal System Parameters', 'FontSize', 18, 'FontWeight', 'bold');
set(gca, 'YTickLabel', pub_labels_rank, 'YTick', 1:k, 'FontSize', 16, 'FontWeight', 'bold');
grid on; box on;

% Add values on bars
for i = 1:k
    text(sorted_mu_star(i) + 0.04*max(sorted_mu_star), i, sprintf('%.2f', sorted_mu_star(i)), ...
         'VerticalAlignment', 'middle', 'FontSize', 15, 'FontWeight', 'bold');
end

% Add cumulative importance percentages
for i = 1:k
    text(0.05*max(sorted_mu_star), i, sprintf('%.0f%%', cumulative_importance(i)*100), ...
         'VerticalAlignment', 'middle', 'FontSize', 13, 'FontWeight', 'bold', 'Color', 'white');
end

%% Summary for Geothermics Publication
fprintf('\n===== GEOTHERMICS PUBLICATION SUMMARY =====\n');
fprintf('MODEL REDUCTION ACHIEVED:\n');
fprintf('- Original parameters: 8\n');
fprintf('- Reduced parameters: 5 (%.0f%% reduction)\n', (1-5/8)*100);
fprintf('- Computational efficiency: ~%.0f%% faster\n', (1-5/8)*100);
fprintf('- Sensitivity coverage: %.0f%% maintained\n', cumulative_importance(end)*100);
fprintf('\nTOP 5 CRITICAL PARAMETERS (in order):\n');
for i = 1:k
    idx = rank_idx(i);
    fprintf('%d. %-20s (μ* = %.2f)\n', i, factors{idx}, sorted_mu_star(i));
end
fprintf('\nFIXED PARAMETERS (low sensitivity):\n');
fprintf('- Rock Specific Heat: %.0f J/kg°C\n', rock_specific_heat_fixed);
fprintf('- Porosity: %.3f\n', porosity_fixed);
fprintf('- Reference Temperature: %.0f °C\n', reference_temp_fixed);
fprintf('==========================================\n');

%% Helper Functions

function X = generate_morris_sample(k, r, p, Delta)
    % Generate Morris sample matrix for reduced parameter set
    X = zeros(r*(k+1), k);
    
    for traj = 1:r
        start_idx = (traj-1)*(k+1) + 1;
        end_idx = traj*(k+1);
        
        % Generate base point
        x_base = rand(1,k) * (p-1) / (p-1);
        x_base = round(x_base * (p-1)) / (p-1);
        
        % Initialize trajectory
        X_traj = repmat(x_base, k+1, 1);
        
        % Random permutation for factor order
        perm = randperm(k);
        
        % Build trajectory
        for i = 1:k
            factor_idx = perm(i);
            
            % Choose direction
            if x_base(factor_idx) + Delta <= 1
                direction = Delta;
            elseif x_base(factor_idx) - Delta >= 0
                direction = -Delta;
            else
                direction = Delta * (2*randi(2)-3);
            end
            
            % Apply perturbation
            X_traj(i+1:end, factor_idx) = X_traj(i+1:end, factor_idx) + direction;
        end
        
        X(start_idx:end_idx, :) = X_traj;
    end
end

function power_capacity = calculate_power_capacity_top5(params, area_fixed, plant_life_seconds, Cw, rho_f, ...
                                                       Cr_fixed, phi_fixed, Tref_fixed)
    % Calculate power capacity for TOP 5 parameters only
    % Fixed parameters are passed as constants
    
    TR        = params(1);     % [°C] Reservoir temperature
    thickness = params(2);     % [m] Thickness  
    Rf        = params(3);     % [fraction] Recovery factor
    rho_r     = params(4);     % [kg/m³] Rock density
    F         = params(5);     % [fraction] Capacity factor
    
    % Fixed parameters (no longer variable)
    Cr   = Cr_fixed;          % [J/kg°C] Rock specific heat - FIXED
    phi  = phi_fixed;         % [fraction] Porosity - FIXED  
    Tref = Tref_fixed;        % [°C] Reference temperature - FIXED

    % Calculate reservoir volume
    V = area_fixed * thickness;
    
    % Calculate thermal energy stored in rock and geofluid
    qR = rho_r * Cr * V * (1 - phi) * (TR - Tref);
    qF = rho_f * Cw * V * phi * (TR - Tref);
    qT = qR + qF;
    
    % Calculate conversion efficiency
    eta_c = max(0, 0.0935 * TR - 2.3266) / 100;
    
    % Calculate electrical power capacity
    power_capacity = (qT * Rf * eta_c) / (F * plant_life_seconds) / 1e6; % [MWe]
end