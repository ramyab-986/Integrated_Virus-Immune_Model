close all; clear; clc;

%% Load data
load('session.mat', 'SS');
load('param.mat', 'param');

iter_prm = size(SS, 3);   % number of parameter sets
iter_ini = size(SS, 1);   % number of initial condition sets

%% Variable indices of interest
var_list   = [1, 70];                  % 1 = VT, 70 = ISGav
var_names  = {'VT', 'ISGav'};          % corresponding variable names

%% Preallocate matrices
VT     = zeros(iter_prm, iter_ini);
ISGav  = zeros(iter_prm, iter_ini);

%% Extract values (log10 transformed)
for ind_prm = 1:iter_prm
    for ind_ini = 1:iter_ini
        VT(ind_prm, ind_ini)    = log10(SS(ind_ini, var_list(1), ind_prm));
        ISGav(ind_prm, ind_ini) = log10(SS(ind_ini, var_list(2), ind_prm));
    end
end

%% Save results
save('VT_ISGav.mat', 'VT', 'ISGav');

%% Reshape data into vectors for visualization
viral_titre = VT(:);
isg_values  = ISGav(:);

%% Plot histogram of viral titre

% Fit a Gaussian Mixture Model with 2 components
gmm_model = fitgmdist(viral_titre, 2);

% Cluster the data based on the GMM
idx = cluster(gmm_model, viral_titre);

% Get proportions of each component
component_counts = histcounts(idx, 2); % Since there are 2 components
total_counts = sum(component_counts);
percentages = 100 * component_counts / total_counts;

fprintf('Percentage in High titre: %.2f%%\n', percentages(1));
fprintf('Percentage in Low titre: %.2f%%\n', percentages(2));

figure;
histogram(viral_titre, 'Normalization', 'pdf', 'EdgeColor', 'none');
set(gca, 'FontSize', 18, 'LineWidth', 1.5, 'box', 'off', 'TickDir', 'both')
xlabel('Viral titre (log_{10})');
ylabel('Probability Density');
title('Distribution of Viral load');
hold on;

% Generate x values for plotting the GMM pdfs
x_vals = linspace(min(viral_titre), max(viral_titre), 1000)';
pdf_vals = pdf(gmm_model, x_vals);
plot(x_vals, pdf_vals, 'r-', 'LineWidth', 2);
legend('Data Histogram', 'GMM Fit', 'Location', 'northwest', 'box', 'off');

%% Correlation plotting

close all; clear; clc;

load("bounds.mat", "PARAMETER")
load('session.mat', 'SS');
load('param.mat', 'param');
load('VT_ISGav.mat', 'VT', 'ISGav')

iter_prm = size(SS, 3);
iter_ini = size(SS, 1);

viral_titre = VT(:);
isg_values  = ISGav(:);

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
    "k_{en, V}", "k_{f, V}", "\mu_{V_{I}}", "k_{a, V}", "k_{e, V}", "k_{r, V}", "k_{c, V}", "k_{t, V}", "\tau", "k_{l, V}", ...
    "N_C", "\mu_r", "\mu_p", "\mu_{V_{E}}", "nSP", "k_s", "k71", "k79", "k76", ...
    "\mu_IRF7", "k72", "k74", "k73", "k75", "tau_5", "k78", "k77", "k70", ...
    "k_{transISGn}", "k69", "degRecBySOCS", "degARCBySOCS", "kinhBySOCS"
];

for ind = 1:129
    C(ind, 1) = corr(PARAMETER(:, ind), viral_titre);
end 

valid_idx = ~isnan(C);
C = C(valid_idx);
param_names = param_names(valid_idx);

[Sorted_C, Sorted_indices] = sort(C);

Names = param_names(Sorted_indices);

% bar(Names, Sorted_C)

% Select top 10 and bottom 10
top10_vals = Sorted_C(end-9:end);
top10_names = Names(end-9:end);

bottom10_vals = Sorted_C(1:10);
bottom10_names = Names(1:10);

% Combine for plotting
final_vals = [bottom10_vals; top10_vals];
final_names = [bottom10_names top10_names];

% Bar plot
figure;
bar(final_names, final_vals)
set(gca, 'XTickLabelRotation', 90, 'LineWidth', 1.5, 'FontSize', 14)
ylabel('Correlation with Viral Load')
