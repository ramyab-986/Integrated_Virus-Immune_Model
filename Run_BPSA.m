%% Bifurcation Point Sensitivity Analysis
% This script performs three analyses:
% 1. Varying each parameter in 'param' (except param 1 and 2; volume ratios- Vc2n and Vn2c are not being changed).
% 2. Computing sensitivity of viral antagonism using a five-point central finite difference.
% 3. Plotting top parameters 
%
% This is extended to immune activity values by appropriate modification. 

clear; close all; clc;

%% Load Data
load('param.mat');
load('SteadyState_120h.mat', 'Tss', 'Yss');

y0_base = Yss(end,:);   % Baseline steady-state values

%% Parameter Names
param_names = [
    "Vc2n", "Vn2c", "b_{IRF3}", "b_{MAVS}", "b_{RIGI}", ...
    "b_{KINASE}", "k1", "k10", "k11", "k12", "k13", "k14", "k15", "k16", "k17", "k18", ...
    "k19", "k2", "k20", "k21", "k22", "k23", "k27", "k29", "k3", "k31", "k32", "k34", ...
    "k35", "k36", "k37", "k38", "k39", "k4", "k40", "k41", "k42", "k43", "k44", "k45", "k46", "k47", "k48", ...
    "k49", "k5", "k50", "k51", "k52", "k53", "k54", "k56", "k57", "k58", ...
    "k59", "k6", "k60", "k61", "k62", "k63", "k64", "k65", "k66", "k67", ...
    "k7", "k8", "k9", "k_{IFN}", "k_{IKK}", "k_{IKKe-TBK1}", "k_{IRF3-IKKe-TBK1}", ...
    "k_{MAVS}", "k_{RIGI}", "k_TFBS_IFNa", "k_TFBS_IFNb", ...
    "k_TFBS_IFNl", "k_act", "k_deph", "k_expr_IkBa", "k_inh_p65", ...
    "k_mRNA_IFNb", "k_mRNA_IFNl", ...
    "k_{mISG}", "k_ISG_mRNA", "k_{rigi synt}", "k_{t, ISG RNA}", "k_trans_IFNl", ...
    "k_transp_NFkB", "mu_IFN", "mu_IFNl", "mu_IkBa", "mu_mRNA_IFNa", ...
    "\mu_mRNA_IFNb", "\mu_mRNA_IFNl", "\mu_{ISG RNA}", "\mu_{ISG}", "\mu_{rigi}", ...
    "k_{en, V}", "k_{f, V}", "\mu_{V_{I}}", "k_{a, V}", "k_{e, V}", "k_{r, V}", "k_{c, V}", "k_{t, V}", "\\tau", "k_{l, V}", ...
    "N_C", "\mu_r", "\mu_p", "\mu_{V_{E}}", "nSP", "k_s", "k71", "k79", "k76", ...
    "\mu_IRF7", "k72", "k74", "k73", "k75", "tau_5", "k78", "k77", "k70", ...
    "k_{transISGn}", "k69", "degRecBySOCS", "degARCBySOCS", "kinhBySOCS"
];

%% Settings
scaling_factors = [0.9, 0.95, 1, 1.05, 1.1, 1.25];
VC_logspace = linspace(log10(0.03), log10(10), 501);
VC_values = 10.^VC_logspace;  % Sweep for V_{I, N}

% Simulation times
t_start = 0;
t_end   = 216 * 60;  % 216 hours in minutes

% Global-like variables
I_a = 1.01388;
I_n = 1;

%% =====================================
% Parameter Sensitivity Analysis (excluding param 1 and 2)
%% =====================================

y0 = Yss(end,:);
y0(2) = 100;   % Fix MOI = 100

critical_VC_data = struct();

disp('--- Starting Parameter sensitivity analysis ---')
tic
for p = 1:length(param)
    % Skip param 1 and 2
    if p == 1 || p == 2
        continue;
    end

    original_param_value = param(p);

    for i = 1:length(scaling_factors)
        param(p) = original_param_value * scaling_factors(i);
        variable_values = zeros(1, length(VC_values));

        % Parallel sweep across VC values
        parfor j = 1:length(VC_values)
            VC = VC_values(j); 
            tspan = linspace(t_start, t_end);
            [~, Y] = ode23s(@(t, y) ODEs(t, y, param, I_n, I_a, VC), tspan, y0);
            variable_values(j) = Y(end, 1);
        end

        % Compute differential wrt log(VC)
        differential = diff(log10(variable_values)) ./ diff(log10(VC_values));
        [~, max_idx] = max(differential);
        critical_VC = VC_values(max_idx + 1);

        % Store results with descriptive parameter name
        param_name = param_names(p);
        critical_VC_data.(matlab.lang.makeValidName(param_name)).(['a_', num2str(scaling_factors(i)*100)]) = critical_VC;
    end

    % Restore parameter
    param(p) = original_param_value;
end
toc

%% =====================================
% 3. Central Difference Sensitivity Calculation
%% =====================================

central_diff_results = struct();

all_param_names = param_names;

for p = 1:length(all_param_names)
    param_name = all_param_names{p};
    safe_name = matlab.lang.makeValidName(param_name);

    if  isfield(critical_VC_data, safe_name)
        critical_VC_0_9  = critical_VC_data.(safe_name).a_90;
        critical_VC_0_95 = critical_VC_data.(safe_name).a_95;
        critical_VC_1_05 = critical_VC_data.(safe_name).a_105;
        critical_VC_1_1  = critical_VC_data.(safe_name).a_110;
    else
        continue; % Skip param_1 and param_2 (not varied)
    end

    % Five-point central difference formula
    central_difference = (critical_VC_0_9 - (8 * critical_VC_0_95) + (8 * critical_VC_1_05) - critical_VC_1_1) / (12 * 0.05);
    central_diff_results.(safe_name) = central_difference;
end

%% Normalization
param_fields = fieldnames(central_diff_results);
sensitivity_values = zeros(1, length(param_fields));
for p = 1:length(param_fields)
    sensitivity_values(p) = central_diff_results.(param_fields{p}) / critical_VC_data.(param_fields{p}).a_100;
end

%% Save All Results
save('Sensitivity.mat', 'central_diff_results');
save('critical_VC_data_combined.mat', 'critical_VC_data');
save('sensitivity_values_normalized_VC_HCV', "sensitivity_values", "param_fields");

%% Plotting top parameters
clear; clc; close all; 

param = [
    "b_{IRF3}", "b_{MAVS}", "b_{RIGI}", ...
    "b_{KINASE}", "k1", "k10", "k11", "k12", "k13", "k14", "k15", "k16", "k17", "k18", ...
    "k19", "k2", "k20", "k21", "k22", "k23", "k27", "k29", "k3", "k31", "k32", "k34", ...
    "k35", "k36", "k37", "k38", "k39", "k4", "k40", "k41", "k42", "k43", "k44", "k45", "k46", "k47", "k48", ...
    "k49", "k5", "k50", "k51", "k52", "k53", "k54", "k56", "k57", "k58", ...
    "k59", "k6", "k60", "k61", "k62", "k63", "k64", "k65", "k66", "k67", ...
    "k7", "k8", "k9", "k_{IFN}", "k_{IKK}", "k_{IKKe-TBK1}", "k_{IRF3-IKKe-TBK1}", ...
    "k_{MAVS}", "k_{RIGI}", "k_TFBS_IFNa", "k_TFBS_IFNb", ...
    "k_TFBS_IFNl", "k_act", "k_deph", "k_expr_IkBa", "k_inh_p65", ...
    "k_mRNA_IFNb", "k_mRNA_IFNl", ...
    "k_{mISG}", "k_ISG_mRNA", "k_{rigi synt}", "k_{t, ISG RNA}", "k_trans_IFNl", ...
    "k_transp_NFkB", "mu_IFN", "mu_IFNl", "mu_IkBa", "mu_mRNA_IFNa", ...
    "\mu_mRNA_IFNb", "\mu_mRNA_IFNl", "\mu_{ISG RNA}", "\mu_{ISG}", "\mu_{rigi}", ...
    "k_{en, V}", "k_{f, V}", "\mu_{V_{I}}", "k_{a, V}", "k_{e, V}", "k_{r, V}", "k_{c, V}", "k_{t, V}", "\\tau", "k_{l, V}", ...
    "N_C", "\mu_r", "\mu_p", "\mu_{V_{E}}", "nSP", "k_s", "k71", "k79", "k76", ...
    "\mu_IRF7", "k72", "k74", "k73", "k75", "tau_5", "k78", "k77", "k70", ...
    "k_{transISGn}", "k69", "degRecBySOCS", "degARCBySOCS", "kinhBySOCS"
];

load('sensitivity_values_normalized_vc_HCV.mat')
sign_vec = 2*(sensitivity_values>0)-1;
[SI_vec, sort_idx] = sort(abs(sensitivity_values), 'descend');
rearr_sign = sign_vec(sort_idx);

new = SI_vec>0.5; 
new_SI = SI_vec(new);
new_sign = rearr_sign(new);
param_sorted = param(sort_idx); 
param_filtered = param_sorted(new);

colors = zeros(length(new_sign), 3); 
colors(new_sign == 1, :) = repmat([0.9, 0.9, 0.9], sum(new_sign == 1), 1); % Red for positive
colors(new_sign == -1, :) = repmat([0, 0, 0], sum(new_sign == -1), 1); % Black for negative

figure;
b = barh(new_SI, 'FaceColor','flat');
b.CData = colors;
set(gca,'FontSize',14,'LineWidth',1.5,'TickDir','both', 'TickLabelInterpreter','tex','XAxisLocation','top', 'YTick', 1:length(param_filtered), 'YTickLabel', param_filtered)%, 'YTickLabelRotation',90);
ylabel('Parameters');
xlabel('Sensitivity Index');
title('HCV - C_V');
box("off")
saveas(gcf, 'Sensitivity_sorted.fig');


% End of script
