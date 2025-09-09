%% Critical Thresholds Analysis under Perturbed Parameter Conditions

%% ========================================================================
% Script to compute critical immune activation values I_{V, N} (for simplicity, we name it as Ia) for 
% different viral antagonism values V_{I, N} (for simplicity, we name it as VC) and parameter perturbations in the 
% ODEs model.
%
% Description:
% - Loads parameter set and steady-state values.
% - Runs ODE simulations across a range of VC.
% - Identifies the corresponding critical I_a value.
% - Repeat analysis for parameter perturbations (0.5× and 2×) for the following parameters:
%         - k_{t, ISG RNA}, mu_{ISG RNA}, k_{r, V}, k_{c, V}, k_{t, V}
% - Saves results to Excel files for each configuration.
% ========================================================================
clear; close all; clc;

% ======================
% Load parameters and steady-state
% ======================
pr = load('param.mat');
load('SteadyState_120h.mat', 'Tss', 'Yss');

y0 = Yss(end,:);
y0(2) = 100;           % MOI = 10

tend = 216 * 60;       % total simulation time (minutes)
tspan = linspace(0, tend);

param_indices = [85, 94, 102, 103, 104]; 
param_names = {'k_{t, ISG RNA}', 'mu_{ISG RNA}', ...
               'k_{r, V}', 'k_{c, V}', 'k_{t, V}'}; 

c = linspace(log10(0.001), log10(5000), 501);
Ia_values = 10.^c;

vc_values = 10.^linspace(-2, 1, 21);

tic
% ======================
% Wild Type
% ======================
disp('Running WT ...');
critical_values_wt = zeros(length(vc_values), 2);

parfor i = 1:length(vc_values)
    vc_local = vc_values(i);
    variable_values = [];
    I_n = 1;

    for j = 1:size(c, 2)
        I_a = Ia_values(j); 
        [T, Y] = ode23s(@(t, y) ODEs(t, y, pr.param, ...
            I_n, I_a, vc_local), tspan, y0);

        variable_values = [variable_values, Y(end, 1)];
    end

    differential = diff(log10(variable_values)) ./ diff(log10(Ia_values));
    [~, max_idx] = max(-differential);
    critical_I_a = Ia_values(max_idx + 1);

    critical_values_wt(i,:) = [vc_local, critical_I_a];
end

wt_excel_filename = fullfile(pwd, ['CriticalValues_WT_', '.xlsx']);
writecell({'VC', 'I_A'}, wt_excel_filename, 'Sheet', 'CriticalValues', 'Range', 'A1');
writematrix(critical_values_wt, wt_excel_filename, 'Sheet', 'CriticalValues', 'Range', 'A2');


% ======================
% Parameter perturbations (0.5× and 2×)
% ======================
for param_idx = 1:length(param_indices)
    param_name = param_names{param_idx};
    disp(['Running param: ', param_name]);

    critical_values_half = zeros(length(vc_values), 2);
    critical_values_double = zeros(length(vc_values), 2);

    parfor i = 1:length(vc_values)
        disp(i)
        vc_local = vc_values(i);
        I_n = 1;

        % ---- Half param ----
        variable_values_half = [];
        param_original = pr.param(param_indices(param_idx));
        pr_half = pr; 
        pr_half.param(param_indices(param_idx)) = param_original * 0.5;

        for j = 1:size(c, 2)
            I_a = Ia_values(j);
            [T, Y] = ode23s(@(t, y) ODEs(t, y, pr_half.param, ...
                I_n, I_a, vc_local), tspan, y0);
            variable_values_half = [variable_values_half, Y(end, 1)];
        end

        diff_half = diff(log10(variable_values_half)) ./ diff(log10(Ia_values));
        [~, max_idx_half] = max(-diff_half);
        critical_Ia_half = Ia_values(max_idx_half + 1);
        critical_values_half(i,:) = [vc_local, critical_Ia_half];

        % ---- Double param ----
        variable_values_double = [];
        pr_double = pr; 
        pr_double.param(param_indices(param_idx)) = param_original * 2.0;

        for j = 1:size(c, 2)
            I_a = Ia_values(j);
            [T, Y] = ode23s(@(t, y) ODEs(t, y, pr_double.param, ...
                I_n, I_a, vc_local), tspan, y0);
            variable_values_double = [variable_values_double, Y(end, 1)];
        end

        diff_double = diff(log10(variable_values_double)) ./ diff(log10(Ia_values));
        [~, max_idx_double] = max(-diff_double);
        critical_Ia_double = Ia_values(max_idx_double + 1);
        critical_values_double(i,:) = [vc_local, critical_Ia_double];
    end

    % Save results
    half_excel_filename = fullfile(pwd, ...
        ['CriticalValues_half_', param_name, '_', '.xlsx']);
    writecell({'VC', 'I_A'}, half_excel_filename, 'Sheet', 'CriticalValues', 'Range', 'A1');
    writematrix(critical_values_half, half_excel_filename, 'Sheet', 'CriticalValues', 'Range', 'A2');

    double_excel_filename = fullfile(pwd, ...
        ['CriticalValues_double_', param_name, '_', '.xlsx']);
    writecell({'VC', 'I_A'}, double_excel_filename, 'Sheet', 'CriticalValues', 'Range', 'A1');
    writematrix(critical_values_double, double_excel_filename, 'Sheet', 'CriticalValues', 'Range', 'A2');
end

toc


%% ===============
% plotting
% ================

% Define parameter names and corresponding filenames
param_names_half = {'half_k_{t, ISG RNA}', 'half_mu_{ISG RNA}', 'half_k_{r, V}', 'half_k_{c, V}', 'half_k_{t, V}'};
param_names_double = {'double_k_{t, ISG RNA}', 'double_mu_{ISG RNA}', ...
                      'double_k_{r, V}', 'double_k_{c, V}', 'double_k_{t, V}'};
param_names = {'k_{t, ISG RNA}', '\mu_{ISG RNA}', 'k_{r, V}', 'k_{c, V}', 'k_{t, V}'};

file_prefix = 'CriticalValues_';
Colors =  [0 0.45 0.74; ...
          0.85 0.33 0.10];

for i = 1:length(param_names_half)
    param_name = param_names_half{i};
    param_name1 = param_names_double{i};
    
    % Construct filenames
    wt_file = fullfile(pwd, ['CriticalValues_WT_', '.xlsx']);
    half_file = fullfile(pwd, [file_prefix, param_name, '_', '.xlsx']);
    double_file = fullfile(pwd, [file_prefix, param_name1, '_', '.xlsx']);
    
    % Load data from Excel
    wt_data = readmatrix(wt_file, 'Sheet', 'CriticalValues', 'Range', 'A2:B100');
    half_data = readmatrix(half_file, 'Sheet', 'CriticalValues', 'Range', 'A2:B100');
    double_data = readmatrix(double_file, 'Sheet', 'CriticalValues', 'Range', 'A2:B100');
    
    % Extract values
    VC_wt = wt_data(:,1);      I_A_wt = wt_data(:,2);
    VC_half = half_data(:,1);  I_A_half = half_data(:,2);
    VC_double = double_data(:,1);  I_A_double = double_data(:,2);
    
    % Plot
    figure(i);
    hold on;
    plot(VC_wt, I_A_wt, '-', 'DisplayName', 'WT', 'LineWidth', 2, 'Color', Colors(1,:));
    plot(VC_half, I_A_half, '--', 'DisplayName', 'Half Param', 'LineWidth', 2, 'Color', Colors(1,:));
    plot(VC_double, I_A_double, 's-', 'DisplayName', 'Double Param', 'LineWidth', 2, 'Color', Colors(1,:));
    
    % Formatting
    set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 20, 'LineWidth', 2, ...
             'YMinorTick', 'off', 'XMinorTick', 'off', 'TickDir', 'both');
    xlabel('VC');
    ylabel('I_A');
    title(param_names{i}, 'Interpreter', 'tex', 'FontSize', 22);
    legend('Location', 'northeastoutside', 'FontSize', 16, 'Box','off');
    grid off;
end
