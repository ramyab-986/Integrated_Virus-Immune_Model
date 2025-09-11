%% Description
% This script simulates and analyzes the effects of sequential IFN (Interferon) 
% dosing on an JAK-STAT signaling model.
%
% Workflow:
% 1. Load steady-state conditions and model parameters.
% 2. Apply a first IFN dose (f1) at time 0 and simulate for 24h.
% 3. Apply a second IFN dose (fr) at 24h and simulate for another 48h.
%
% Analysis:
% - Extract time courses of selected variables involved in JAK-STAT signaling.
% - Compute chi-squared metrics to quantify differences between two IFN dosing regimens.

%% IFN exposure simulation
clear; close all; clc;

load('param.mat');
load('SteadyState_120h.mat', 'Tss', 'Yss');

y0_ss = Yss(end,:);  % Steady-state initial condition

I_n = 1; VC = 0; I_a = 0;

f1 = [0.1, 1, 10, 100];
fr = [0.1, 1, 10, 100];

sst = 24*1*60;        %1st IFN dose exposure time
tspan1 = linspace(0, sst);

simt = 24*2*60;     %2nd IFN dose exposure time
tspan2 = linspace(0, simt);

for i = 1:length(f1)
    
    y0 = y0_ss;
    y0(39) = f1(i);  %IFN input

    [T1, Y1] = ode23s(@(t, y) ODEs(t, y, param, I_n, I_a, VC), tspan1, y0);
    save(['IFN_24h_', num2str(f1(i)), '.mat'], 'T1', 'Y1');

    for j = 1:length(fr)
        initial_values = Y1(end,:);
        initial_values(39) = initial_values(39) + fr(j);

        [T2, Y2] = ode23s(@(t, y) ODEs(t, y, param, I_n, I_a, VC), tspan2, initial_values);
        save(['IFN_', num2str(f1(i)), '_IFNrest', num2str(fr(j)), '_48h.mat'], 'T2', 'Y2');

        combined_T = [T1; T2 + T1(end)];
        combined_Y = [Y1; Y2];
        save(['combined_IFN_', num2str(f1(i)), '_IFNrest', num2str(fr(j)), '_48h.mat'], 'combined_T', 'combined_Y');
    end
end

%%
% Parameters
f1 = [1, 100];   % First dose values
f2 = [1, 100];   % Second dose values
lw = 2;          % Line width

% Define markers and colors
symbols = {'o','s','^','d','v','p','h','x','*','+','.'};
colors = [0.93, 0.69, 0.13;
          0.85, 0.33, 0.10];

% Variable names
var_names = {'ExtVirus','VirusInit','IntVirus','R_{cyt}','(+)RNA_{CM}','SP','NSP','RC_{CM}','dsRNA','RIGI','aRIGI','MAVS','aMAVS',...
    'IKKe','pIKKe','TBK1','pTBK1','IRF3','pIRF3','IKK','aIKK','NFkBIkBac','pNFkBn','NFkBn','NFkBc','IkBac','IRF7','pIRF7','IFNbmRNA',...
    'IFNamRNA','IFNlmRNA','IFN_c','IFNl_c','JAK','RJC','STAT1c','CP','ISGn','IFNex','STAT2c','TYK','RTC','ARC','Rec1','Rec2',...
    'IFNARd','IRF9_c','ARC-STAT2_c','ARC-STAT12_c','STAT2-IRF9_c','ISGF3_c','PSC_c','ISGF3-CP','PSC-CP','NP','STAT1_n','STAT2_n','PIAS','PSC_n',...
    'IRF9_n','ISGF3_n','PSC-NP','B_u','B_o-NP','B_o','ISGF3-PIAS','STAT2-IRF9_n','ISGF3n-NP','ISGavmRNA','ISGav','ISGnmRNA_n','IRF9mRNA_n',...
    'IRF7mRNA','ISGnmRNA_c','IRF9mRNA_c'};

% Variables of interest
v = [69, 74];  %ISGavmRNA, IRF9mRNA_c

% Loop over variables
for k = 1:length(v)
    var_idx = v(k);

    % Loop over first dose values
    for i = 1:length(f1)
        global_ymax = 0; % reset per plot

        % New figure
        figure; hold on;
        title(sprintf('%s, D_1 = %dnM', var_names{var_idx}, f1(i)), 'FontSize', 20);

        % Load first dataset
        filename1 = sprintf('IFN_24h_%d.mat', f1(i));
        data1 = load(filename1);

        % Update ymax
        global_ymax = max(global_ymax, max(data1.Y1(:, var_idx)));

        % Plot first dose
        plot(data1.T1 / 60, data1.Y1(:, var_idx), ...
            'LineWidth', lw, ...
            'DisplayName', sprintf('D_1 = %d only', f1(i)));

        % Loop over second dose values
        for j = 1:length(f2)
            filename2 = sprintf('IFN_%d_IFNrest%d_48h.mat', f1(i), f2(j));
            data2 = load(filename2);

            % Update ymax
            global_ymax = max(global_ymax, max(data2.Y2(:, var_idx)));

            % Plot second dose
            plot((data2.T2 / 60) + 24, data2.Y2(:, var_idx), ...
                'LineWidth', lw, ...
                'Color', colors(j, :), ...
                'Marker', symbols{j}, ...
                'MarkerIndices', 1:10:length(data2.T2), ...
                'MarkerSize', 8, ...
                'DisplayName', sprintf('D_1 = %d, D_2 = %d', f1(i), f2(j)));
        end

        % Mark 24h boundary
        xline(24, 'k--', 'LineWidth', 2);

        % Formatting
        set(gca, 'FontSize', 18, 'LineWidth', 1.5, 'TickDir', 'both');
        xlabel('Time (h)', 'FontSize', 18);
        ylabel(var_names{var_idx}, 'FontSize', 18);
        xlim([-2, 74]);
        ylim([0, global_ymax]);

        % Legend
        legend('Location', 'northeastoutside', 'AutoUpdate', 'off');
        legend boxoff;
    end
end

%% Chi-sqaured Analysis and Plotting

clear; close all; clc;

f1 = [1, 100];  % 1st IFN dose
f2_r = [1, 100];  %2nd IFN dose

var_names = {'ExtVirus', 'VirusInit', 'IntVirus', 'R_{cyt}', '(+)RNA_{CM}', 'SP', 'NSP', 'RC_CM', 'dsRNA', 'RIGI','aRIGI','MAVS','aMAVS',...
    'IKKe','pIKKe','TBK1','pTBK1', 'IRF3','pIRF3','IKK','aIKK','NFkBIkBac','pNFkBn','NFkBn','NFkBc', 'IkBac', 'IRF7', 'pIRF7', 'IFNbmRNA',...
    'IFNamRNA','IFNlmRNA', 'IFN_c', 'IFNl_c', 'JAK','RJC', 'STAT1c','CP', 'ISGn','IFNex','STAT2c','TYK','RTC','ARC', 'Rec1','Rec2',...
    'IFNARd','IRF9_c','ARC-STAT2_c', 'ARC-STAT12_c','STAT2-IRF9_c','ISGF3_c', 'PSC_c','ISGF3-CP','PSC-CP','NP','STAT1_n','STAT2_n','PIAS','PSC_n',...
    'IRF9_n','ISGF3_n','PSC-NP','B_u','B_o-NP','B_o','ISGF3-PIAS','STAT2-IRF9_n','ISGF3n-NP', 'ISGavmRNA','ISGav', 'ISGnmRNA_n', 'IRF9mRNA_n',...
    'IRF7mRNA', 'ISGnmRNA_c', 'IRF9mRNA_c'};

ind = [10, 27, 34:38, 40:75]; 
rowLabels = var_names(ind);

data = struct([]);

for i = 1:length(f1)
    filename = ['combined_IFN_', num2str(f1(i)), '_IFNrest', num2str(100), '_48h.mat'];
    if isfile(filename)
        tempData = load(filename);
        if isempty(data)
            data = tempData;
        else
            data(i) = tempData;
        end
    end
end

arr = zeros(200, 43, length(f1));
for i = 1:43
    for j = 1:length(f1)
        arr(:, i, j) = reshape([data(j).combined_Y(:, ind(i))], [], 1);
    end
end

time = data(1).combined_T / 60;  
numVariables = size(arr, 2);
numCases = length(f1);

% Chi-squared computation with log10 safety
epsilon = 1e-6;
chiSquaredValues = zeros(numVariables, nchoosek(numCases, 2));
caseIdx = 1;
for i = 1:numCases-1
    for j = i+1:numCases
        for k = 1:numVariables
            vals1 = arr(:, k, i);
            vals2 = arr(:, k, j);

            % Replace 0s with epsilon where necessary
            bothZeroIdx = vals1 == 0 & vals2 == 0;
            vals1(bothZeroIdx) = epsilon;
            vals2(bothZeroIdx) = epsilon;

            vals1(vals1 == 0) = epsilon;
            vals2(vals2 == 0) = epsilon;

            obs1 = log10(vals1);
            obs2 = log10(vals2);

            chiSquaredValues(k, caseIdx) = (sum((obs1 - obs2).^2))/200;
        end
        caseIdx = caseIdx + 1;
    end
end

avgChiSquaredValues = mean(chiSquaredValues, 2);
[sort_ChiSquaredValues, sort_idx] = sort(avgChiSquaredValues, 'descend');

figure;
bar(sort_ChiSquaredValues);
xlabel('Variable Index');
ylabel('Average \chi^2');
set(gca, 'XTick', 1:1:43, 'XTickLabel', rowLabels(sort_idx), 'XTickLabelRotation', 45);
grid off;

disp('Top 3 variables with highest average chi-squared values:');
disp(rowLabels(sort_idx(1:3)));

% Highlight selected variable time courses
highlight_idx = [3, 9, 10, 14];  

for k = 1:length(highlight_idx)
    idx = highlight_idx(k);
    label = rowLabels(idx);

    figure;
    semilogy(time, arr(:, idx, 1), 'LineWidth', 2); hold on;
    semilogy(time, arr(:, idx, 2), 'LineWidth', 2);
    set(gca, 'FontSize', 18, 'XMinorTick', 'off', 'YMinorTick', 'off', ...
        'XTick', [0, 24, 48, 72], 'LineWidth', 1.5, 'box', 'off', ...
        'TickDir', 'both')
    legend('D_1=1 nM', 'D_1=100 nM', 'Location', 'bestoutside', 'box', 'off');
    xlim([-2 74]);
    xlabel('Time (hours)');
    ylabel('Concentration (nM)');
    title(label);
    hold off;
end

