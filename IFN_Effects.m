%% =============================================================
%  IFN-HCV Simulation
%  Simulates antiviral effects of IFN treatment in three scenarios:
%   1) Pre-infection IFN treatment
%   2) IFN at time of infection (0 h)
%   3) Post-infection IFN treatment
%  Results are computed for multiple viral antagonism (VC here) values and 
%  visualized as heatmaps of final viral loads.
%% =============================================================

clear; close all; clc;

%% ================== Setup ==================
load('param.mat');
load('SteadyState_120h.mat', 'Tss', 'Yss');

IFN_conditions = [0, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30];   % IFN concentrations (nM)
VC_values = [0.25*0.439, 0.5*0.439, 0.439, 2*0.439, 4*0.439];
I_n = 1; I_a = 1.01388; 

% Time arrays (minutes)
timearray_pre  = [-720:24:-24, -16, -8, -4] * 60;   % pre-infection
timearray_zero = 0;                                % at infection
timearray_post = [4, 8, 16, 24, 48, 72] * 60;      % post-infection

timearray_combined = [timearray_pre, timearray_zero, timearray_post] / 60; % hours

%% ================== Loop over VC values ==================
for v = 1:length(VC_values)
    VC = VC_values(v);
    disp(['Running simulations for VC = ', num2str(VC)]);

    %% ----- Case 1: Pre-infection IFN -----
    ResultsPre = struct();
    for i = 1:length(IFN_conditions)
        for k = 1:length(timearray_pre)
            % Initial conditions
            y0 = Yss(end,:);
            y0(39) = IFN_conditions(i); % IFN before infection

            % Simulate until IFN treatment end time
            t_end = timearray_pre(k);
            tspan1 = linspace(0, abs(t_end));
            [~, Y1] = ode23s(@(t, y) ODEs(t, y, param, I_n, I_a, VC), tspan1, y0);

            % Infect with virus
            y01 = Y1(end, :);
            y01(2) = y01(2) + 100; % virus input
            tspan2 = linspace(abs(t_end)+1e-8, abs(t_end)+48*60);
            [~, Y2] = ode23s(@(t, y) ODEs(t, y, param, I_n, I_a, VC), tspan2, y01);

            % Store results
            idx = sub2ind([length(IFN_conditions), length(timearray_pre)], i, k);
            ResultsPre(idx).VC = VC;
            ResultsPre(idx).IFN = IFN_conditions(i);
            ResultsPre(idx).time = timearray_pre(k)/60;
            ResultsPre(idx).ViralLoad = Y2(end,1);
        end
    end
    save(['ResultsPre_VC', num2str(VC), '.mat'], 'ResultsPre');

    %% ----- Case 2: IFN at zero hour -----
    ResultsZero = struct();
    for i = 1:length(IFN_conditions)
        % Initial conditions
        y0 = Yss(end,:);
        y0(2) = 100;                % virus input
        y0(39) = IFN_conditions(i); % IFN at 0h

        % Simulate
        tspan = linspace(0, 48*60);
        [~, Y1] = ode23s(@(t, y) ODEs(t, y, param, I_n, I_a, VC), tspan, y0);

        % Store
        ResultsZero(i).VC = VC;
        ResultsZero(i).IFN = IFN_conditions(i);
        ResultsZero(i).time = 0;
        ResultsZero(i).ViralLoad = Y1(end,1);
    end
    save(['ResultsZero_VC', num2str(VC), '.mat'], 'ResultsZero');

    %% ----- Case 3: Post-infection IFN -----
    ResultsPost = struct();
    for i = 1:length(IFN_conditions)
        for k = 1:length(timearray_post)
            % Initial conditions
            y0 = Yss(end,:);
            y0(2) = 100; % virus input
            y0(39) = 0;  % no IFN initially

            % Simulate until IFN treatment time
            t_end = timearray_post(k);
            tspan1 = linspace(0, t_end);
            [~, Y1] = ode23s(@(t, y) ODEs(t, y, param, I_n, I_a, VC), tspan1, y0);

            % Add IFN
            y01 = Y1(end,:);
            y01(39) = y01(39) + IFN_conditions(i);

            % Continue simulation
            tspan2 = linspace(t_end+1e-8, t_end+48*60);
            [~, Y2] = ode23s(@(t, y) ODEs(t, y, param, I_n, I_a, VC), tspan2, y01);

            % Store
            idx = sub2ind([length(IFN_conditions), length(timearray_post)], i, k);
            ResultsPost(idx).VC = VC;
            ResultsPost(idx).IFN = IFN_conditions(i);
            ResultsPost(idx).time = timearray_post(k)/60;
            ResultsPost(idx).ViralLoad = Y2(end,1);
        end
    end
    save(['ResultsPost_VC', num2str(VC), '.mat'], 'ResultsPost');

    %% ----- Plotting for this VC -----
    % Reconstruct matrices
    viral_load_matrix_pre  = reshape([ResultsPre.ViralLoad],  length(IFN_conditions), length(timearray_pre));
    viral_load_matrix_zero = [ResultsZero.ViralLoad]';
    viral_load_matrix_post = reshape([ResultsPost.ViralLoad], length(IFN_conditions), length(timearray_post));

    viral_load_matrix_combined = [viral_load_matrix_pre, viral_load_matrix_zero, viral_load_matrix_post];

    % Heatmap
    figure('Position', [100, 100, 950, 500]);
    imagesc(timearray_combined, 1:length(IFN_conditions), log10(viral_load_matrix_combined));
    colormap('sky');
    set(gca, 'YTick', 1:length(IFN_conditions), 'YTickLabel', IFN_conditions, ...
        'XTick', -312:17:72, 'XTickLabel', {'-312','-288','-264','-240','-216','-192','-168','-144','-120','-96','-72','-48','-24','-16','-8','-4','0','4','8','16','24','48','72'}, ...
        'FontSize', 20, 'LineWidth', 1.5, 'CLim', [-4.5 4.5]);
    xlim([-312, 72]);

    xlabel('Time (h)');
    ylabel('IFN concentration (nM)');
    c = colorbar;
    set(c, 'Ticks', [-4, -2, 0, 2, 4], 'TickLabels', {'10^{-4}','10^{-2}','10^{0}','10^{2}','10^{4}'}, ...
        'LineWidth', 1.5, 'FontSize', 20);
    c.Label.String = 'Viral Load';
    c.Label.FontSize = 22;

    % Save plots
    saveas(gcf, ['Combined_Heatmap_VC', num2str(VC), '.png']);
    saveas(gcf, ['Combined_Heatmap_VC', num2str(VC), '.fig']);
end
