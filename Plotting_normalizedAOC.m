% =============================================================
%  IFN-HCV nAOC Analysis
%  Loads pre-infection IFN simulation results across multiple 
%  viral antagonism strength (VC here) values. Computes area over curve (AOC) 
%  relative to a viral load threshold during the pre-treatment 
%  window (t ≤ 0), normalizes AOC to a reference VC, and 
%  visualizes the normalized AOC (nAOC) as a heatmap.
% =============================================================

clear; close all; clc;

%% ========================
%  Define time arrays
%  ========================
timearray_pre  = [-720:24:-24, -16, -8, -4] * 60;   % minutes
timearray_combined = timearray_pre / 60; % hours
time_days = timearray_combined / 24;

%% ========================
%  IFN and VC values
%  ========================
IFN_conditions = [0, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30];
VC_values = [0.25*0.439, 0.5*0.439, 0.439, 2*0.439, 4*0.439];

nIFN = length(IFN_conditions);
nTime = length(timearray_combined);
nVC   = length(VC_values);

%% ========================
%  Initialize storage
%  ========================
viral_load_matrix_combined = nan(nIFN, nTime, nVC);
aoc_matrix = nan(nVC, nIFN);

%% ========================
%  Load results from struct files
%  ========================
for v = 1:nVC
    VC = VC_values(v);

    % Load pre/zero/post results
    load(['ResultsPre_VC',  num2str(VC), '.mat'], 'ResultsPre');
    
    % Convert struct arrays into matrices
    viral_load_matrix_pre  = reshape([ResultsPre.ViralLoad],  nIFN, length(timearray_pre));
    
    % Store combined viral load per IFN × VC
    viral_load_matrix_combined(:, :, v) = viral_load_matrix_pre;
end

%% ========================
%  Compute AOC (t ≤ 0) relative to threshold
%  ========================
y_threshold = 3;

for i_IFN = 1:nIFN
    for v = 1:nVC
        viral_load = squeeze(log10(viral_load_matrix_combined(i_IFN, :, v)));
        mask = time_days <= 0;
        diff_from_thresh = y_threshold - viral_load(mask);
        aoc_matrix(v, i_IFN) = trapz(time_days(mask), diff_from_thresh);
    end
end

%% ========================
%  Normalize AOC by reference VC
%  ========================
ref_row = 3; % row corresponding to VC = 0.439
aoc_matrix_n = zeros(size(aoc_matrix(:, 2:end))); % skip IFN=0

for col = 2:nIFN
    aoc_matrix_n(:, col-1) = aoc_matrix(:, col) ./ aoc_matrix(ref_row, col);
end

%% ========================
%  Custom colormap
%  ========================
nums = 15;
map = zeros(nums, 3);
for ind = 1:nums
    t = (ind - 1) / (nums - 1);  
    R = t;
    G = 0.25;
    B = 1 - t * 0.5;
    map(ind, :) = [R, G, B];
end

%% ========================
%  Plot heatmap
%  ========================
figure;
h = heatmap(IFN_conditions(2:end), VC_values, aoc_matrix_n, 'FontSize', 18);

h.Colormap = map;
h.ColorbarVisible = 'on';
h.GridVisible = 'off';
h.CellLabelColor = 'none';

h.XLabel = 'IFN Concentration (nM)';
h.YLabel = 'V_{I,N}/V_{I,N}^0';
h.Title = 'nAOC';
