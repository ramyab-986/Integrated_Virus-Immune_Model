% Description:

%ADD DESCRIPTION

clear; close all; clc;

% Load data
load('param.mat');
load('SteadyState_120h.mat', 'Tss', 'Yss');

% Define parameter and simulation settings
IFN_conditions = [0, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30]; 
timearray = [4, 8, 16, 24:24:720]*60; 

I_n = 1; I_a = 1; VC = 0.439;

t_start = 0;
param_indices = 103;
param_names = {'k_c'};
param_factors = [0.03, 0.1, 0.3, 1, 3, 10, 30];

% Initialize structure to store results
results = struct();

tic
for param_idx = 1:length(param_indices)
    param_name = param_names{param_idx};
    param_original = param(param_indices(param_idx));
    disp(param_original)

    for factor = param_factors
        param(param_indices(param_idx)) = param_original * factor;
        
        for i = 1:length(IFN_conditions)
            y0 = Yss(end,:);
            y0(39) = IFN_conditions(i); % IFN pretreatment

            for k = 1:length(timearray)
                t_end = timearray(k);
                tspan1 = [t_start t_end];
                [T1, Y1] = ode23s(@(t, y) ODEs(t, y, param, I_n, I_a, VC), tspan1, y0);

                y01 = Y1(end, :);
                y01(2) = y01(2) + 100; % Virus input
                t_end1 = t_end + 48*60;
                tspan2 = linspace(t_end + 1e-8, t_end1);
                [T2, Y2] = ode23s(@(t, y) ODEs(t, y, param, I_n, I_a, VC), tspan2, y01);

                % Store results in structure
                factor_str = strrep(num2str(factor), '.', '_'); 
                IFN_str = strrep(num2str(IFN_conditions(i)), '.', '_'); 
                results.(param_name).(['factor_', factor_str]).(['IFN_', IFN_str]).(['time_', num2str(timearray(k)/60), 'h']) = struct('T', [T1; T2], 'Y', [Y1; Y2]);
            end
        end
        param(param_indices(param_idx)) = param_original; % Reset parameter
    end
end
toc

% Save results
save('IFN_Antiviral_kc.mat', 'results');
%% plotting
% 
% load('IFN_Antiviral_kc.mat', 'results');

% Define time points (pre-treatment windows, negative time = before infection)
timearray_pre = [-720:24:-24, -16, -8, -4] * 60;  
timearray_combined = timearray_pre / 60;  

% Define IFN conditions and parameter settings
IFN_conditions = [0, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30];
param_cases = fieldnames(results);    % {'k_c'} here
param_factors = [0.03, 0.1, 0.3, 1, 3, 10, 30];

% Initialize storage
viral_load_data = struct();

for p = 1:length(param_cases)
    param = param_cases{p};

    for i_IFN = 1:length(IFN_conditions)
        IFN = IFN_conditions(i_IFN);
        IFN_key = ['IFN_' strrep(num2str(IFN), '.', '_')]; 

        viral_load_matrix_pre = zeros(length(param_factors), length(timearray_pre));

        for f = 1:length(param_factors)
            factor = param_factors(f);
            factor_key = ['factor_', strrep(num2str(factor), '.', '_')];

            for j = 1:length(timearray_pre)
                time_h = abs(timearray_pre(j)/60);
                time_key = ['time_', num2str(time_h), 'h'];

                sim = results.(param).(factor_key).(IFN_key).(time_key);

                % Take viral load at end of simulation
                viral_load_matrix_pre(f, j) = sim.Y(end, 1);   % viral load
            end
        end

        viral_load_data.(param).(IFN_key) = viral_load_matrix_pre;
    end
end

% Process viral load data -> AOC calculation
y_threshold = 3;  
time_days = timearray_combined / 24;  

aoc_matrix = zeros(length(param_factors), length(IFN_conditions));
for i_IFN = 1:length(IFN_conditions)
    IFN_key = ['IFN_' strrep(num2str(IFN_conditions(i_IFN)), '.', '_')];

    for f_idx = 1:length(param_factors)
        param = param_cases{1};
        viral_load = log10(viral_load_data.(param).(IFN_key)(f_idx, :));

        mask = time_days <= 0; % Only use pre-treatment window
        aoc = trapz(time_days(mask), y_threshold - viral_load(mask));
        aoc_matrix(f_idx, i_IFN) = aoc;
    end
end

% Normalize AOCs by reference row (wildtype factor = 1)
ref_row = find(param_factors == 1);
aoc_matrix_n = zeros(size(aoc_matrix(:, 2:end)));
for r = 1:length(param_factors)
    aoc_matrix_n(r, :) = aoc_matrix(r, 2:end) ./ aoc_matrix(ref_row, 2:end);
end

% Custom red-blue colormap
nums = 15;
map = zeros(nums, 3);
for ind = 1:nums
    t = (ind - 1) / (nums - 1);
    R = t; 
    G = 0.25; 
    B = (1 - 0.5*t);
    map(ind, :) = [R, G, B];
end

% Plot heatmap
figure;
heatmap(IFN_conditions(2:end), param_factors, aoc_matrix_n, ...
     'ColorbarVisible','on', 'GridVisible','off', ...
     'CellLabelColor','none', 'FontSize',18);
colormap(map);
xlabel('IFN (nM)');
ylabel('k_c / k_c^0');
