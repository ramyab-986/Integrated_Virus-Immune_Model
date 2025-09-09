%% ===============================================================
%   Immune Response Dynamics Clustering and Visualization
%
%   This script performs the following steps:
%   1. Loads simulation results (HCV model, VC conditions - critical, 0.5*critical and 2*critical).
%   2. Filters out variables with low variance.
%   3. Normalizes trajectories across conditions (z-score).
%   4. Performs hierarchical clustering (correlation-based).
%   5. Generates multiple visualizations:
%        - Heatmaps with dendrograms.
%        - Cluster mean dynamics across conditions.
%        - Component lists per cluster.
%        - SSD (sum of squared differences) analysis.
%        - Correlation-based distance tables.
%        - Per-component dynamics with cluster markers.
%   6. Saves output tables and figures.
%
% ================================================================

clear; close all; clc;

%% ================== USER OPTIONS ==================
k             = 8;       % number of clusters
vars_range    = 10:75;   % excluding viral life cycle variables
std_thresh    = 0.05;    % relative standard deviation threshold

%% ==== Load data ====
A = load('HCV_48h_VC0.22_Ia1.0138_In1.mat');
B = load('HCV_48h_VC0.439_Ia1.0138_In1.mat');
C = load('HCV_48h_VC0.878_Ia1.0138_In1.mat');

YA = A.Y(:, vars_range);
YB = B.Y(:, vars_range);
YC = C.Y(:, vars_range);

% Variable names (matching vars_range indices)
var_names = { 'RIGI','aRIGI','MAVS','aMAVS', ...
    'IKKe','pIKKe','TBK1','pTBK1', 'IRF3','pIRF3','IKK','aIKK', ...
    'NFkBIkBac','pNFkBn','NFkBn','NFkBc', 'IkBac', 'IRF7', 'pIRF7', ...
    'IFNbmRNA','IFNamRNA','IFNlmRNA', 'IFN_c', 'IFNl_c', 'JAK','RJC', ...
    'STAT1c','CP', 'ISGn','IFNex','STAT2c','TYK','RTC','ARC', 'Rec1','Rec2', ...
    'IFNARd','IRF9_c','ARC-STAT2_c', 'ARC-STAT12_c','STAT2-IRF9_c','ISGF3_c', ...
    'PSC_c','ISGF3-CP','PSC-CP','NP','STAT1_n','STAT2_n','PIAS','PSC_n', ...
    'IRF9_n','ISGF3_n','PSC-NP','B_u','B_o-NP','B_o','ISGF3-PIAS', ...
    'STAT2-IRF9_n','ISGF3n-NP', 'ISGavmRNA','ISGav', 'ISGnmRNA_n', ...
    'IRF9mRNA_n','IRF7mRNA', 'ISGnmRNA_c', 'IRF9mRNA_c'};

%% ==== Filter low-variance + manual removed vars ====
stdA = std(YA); meanA = mean(YA);
stdB = std(YB); meanB = mean(YB);
stdC = std(YC); meanC = mean(YC);

low_std_mask = (stdA./(meanA+eps) < std_thresh) & ...
               (stdB./(meanB+eps) < std_thresh) & ...
               (stdC./(meanC+eps) < std_thresh);


final_mask = ~(low_std_mask);

YA = YA(:, final_mask);
YB = YB(:, final_mask);
YC = YC(:, final_mask);
var_names = var_names(final_mask);

[Tlen, Ngene] = size(YA);
fprintf('Remaining vars after filtering: %d\n', Ngene);

%% ==== Normalize (z-score across conditions) ====
Yall = [YA; YB; YC];
mu    = mean(Yall);
sigma = std(Yall);

zAn = (YA - mu) ./ sigma;
zBn = (YB - mu) ./ sigma;
zCn = (YC - mu) ./ sigma;
z   = [zAn; zBn; zCn];

%% ==== Clustering (concatenated correlation) ====
D_concat_corr    = pdist(z', 'correlation');
Z                = linkage(D_concat_corr, 'average');
T                = cluster(Z, 'maxclust', k);
[~,~,leafOrder]  = dendrogram(Z, 0); 
close(gcf);

%% ==== Visualizations & Analyses ====
plotHeatmaps(Z, leafOrder, var_names, zAn, zBn, zCn, A.T/60, B.T/60, C.T/60);
plotClusterMeans(k, T, zAn, zBn, zCn, A.T);
saveClusterGeneLists(T, var_names);
ssdAnalysis(YA, YC, var_names, T);
distanceTable = computeClusterDistances(k, T, zAn, zBn, zCn, var_names);
writetable(distanceTable, 'cluster_gene_distances.csv');
plotPerClusterDynamics(k, T, zAn, zBn, zCn, var_names, distanceTable.Distance_Mean, A.T);

%% ==== FUNCTIONS ====

function plotHeatmaps(Z, leafOrder, var_names, zAn, zBn, zCn, tA, tB, tC)
    plotScenario('Condition A (Consensus order)', Z, leafOrder, var_names, zAn, tA);
    plotScenario('Condition B (Consensus order)', Z, leafOrder, var_names, zBn, tB);
    plotScenario('Condition C (Consensus order)', Z, leafOrder, var_names, zCn, tC);
end

function plotClusterMeans(k, T, zAn, zBn, zCn, time)
      colors = [
    0.8500    0.3300    0.1000;
    0.9300    0.6900    0.1300;
    0.3000    0.7500    0.9300;
    0.4900    0.1800    0.5600;
    1.0000    0.0700    0.6500;
    0.0000    0.0000    1.0000;
    0.4700    0.6700    0.1900
    0.09      0.81      0.73];

    for c = 1:k
        vars_c = find(T == c);
        mA = mean(zAn(:, vars_c), 2);
        mB = mean(zBn(:, vars_c), 2);
        mC = mean(zCn(:, vars_c), 2);
        figure('Position', [100 100 200 150]);
        lw = 1.5;
        plot(time/60, mA, 'LineWidth', lw, 'LineStyle', ':','Color', colors(c,:)); hold on;
        plot(time/60, mB, 'LineWidth', lw, 'LineStyle', '--', 'Color', colors(c,:));
        plot(time/60, mC, 'LineWidth', lw, 'LineStyle', '-',  'Color', colors(c,:));
        set(gca,'FontSize',16,'LineWidth',1.5,'TickDir','both','XTick',[0 24 48],'Box','off');
        xlim([-1, 50]);
    end
end

function saveClusterGeneLists(T, var_names)
    fid = fopen('vars_with_clusters.txt', 'w');
    fprintf(fid, 'Cluster-wise Gene Components:\n\n');
    for c = 1:max(T)
        vars_c = find(T == c);
        gene_list = var_names(vars_c);
        fprintf('Cluster %d:\n', c);
        fprintf('  %s\n', strjoin(gene_list, ', '));
        fprintf('\n');
        fprintf(fid, 'Cluster %d:\n', c);
        fprintf(fid, '  %s\n\n', strjoin(gene_list, ', '));
    end
    fclose(fid);
end

function ssdAnalysis(YA, YC, var_names, T)
    X1 = log10(YA(1:34,:) + 1e-10);
    X2 = log10(YC(1:34,:) + 1e-10);
    ssd_vals = mean((X1 - X2).^2, 1);
    [sorted_ssd, idx] = sort(ssd_vals, 'descend');
    sorted_vars = var_names(idx);
    
      colors = [
    0.8500    0.3300    0.1000;
    0.9300    0.6900    0.1300;
    0.3000    0.7500    0.9300;
    0.4900    0.1800    0.5600;
    1.0000    0.0700    0.6500;
    0.0000    0.0000    1.0000;
    0.4700    0.6700    0.1900
    0.09      0.81      0.73];

    figure('Position', [100 100 400 250]); hold on;
    for i = 1:length(sorted_ssd)
        cluster_id = T(idx(i));
        edge_color = colors(cluster_id, :);
        bar(i, sorted_ssd(i), 'FaceColor', edge_color, 'EdgeColor', 'none');
    end
    ylabel('\chi^2');
    set(gca,'FontSize',14,'LineWidth',1.5,'TickDir','both','XTick',1:length(sorted_vars),...
        'XTickLabel',sorted_vars,'XTickLabelRotation',90,'Box','off');
end

function distance_table = computeClusterDistances(k, T, zAn, zBn, zCn, var_names)
    distance_table = cell(0,6);
    dist_mean_vec  = nan(1, numel(var_names));
    for c = 1:k
        vars_c = find(T == c);
        mA = mean(zAn(:, vars_c), 2);
        mB = mean(zBn(:, vars_c), 2);
        mC = mean(zCn(:, vars_c), 2);
        for g = 1:numel(vars_c)
            idx_gene = vars_c(g);
            dA = corr(zAn(:, idx_gene), mA, 'rows','complete');
            dB = corr(zBn(:, idx_gene), mB, 'rows','complete');
            dC = corr(zCn(:, idx_gene), mC, 'rows','complete');
            dist_mean = mean([dA,dB,dC]);
            dist_mean_vec(idx_gene) = dist_mean;
            distance_table(end+1,:) = {c, var_names{idx_gene}, dA, dB, dC, dist_mean};
        end
    end
    distance_table = cell2table(distance_table,...
        'VariableNames',{'Cluster','Gene','Distance_A','Distance_B','Distance_C','Distance_Mean'});
end

function plotPerClusterDynamics(k, T, zAn, zBn, zCn, var_names, dist_mean_vec, time)
    colors = [
    0.8500    0.3300    0.1000;
    0.9300    0.6900    0.1300;
    0.3000    0.7500    0.9300;
    0.4900    0.1800    0.5600;
    1.0000    0.0700    0.6500;
    0.0000    0.0000    1.0000;
    0.4700    0.6700    0.1900
    0.09      0.81      0.73];

    time_h = time / 60;
    line_styles = {':','--','-'};
    
    for c =1:k
        vars_c = find(T == c);
        nGc   = numel(vars_c);
        nCols = 4; nRows = ceil(nGc / nCols);
        figure('Name',sprintf('Cluster %d - dynamics',c),'Color','w','Position',[100 100 1200 800]);
              
        for gi = 1:nGc
            idx_gene = vars_c(gi);
            subplot(nRows,nCols,gi);
            
            % Plot trajectories
            plot(time_h, zAn(:, idx_gene),'LineWidth',2,'Color',colors(c,:),'LineStyle',line_styles{1}); hold on;
            plot(time_h, zBn(:, idx_gene),'LineWidth',2,'Color',colors(c,:),'LineStyle',line_styles{2});
            plot(time_h, zCn(:, idx_gene),'LineWidth',2,'Color',colors(c,:),'LineStyle',line_styles{3});
            
            % Add marker encoding distance to centroid
            d = dist_mean_vec(idx_gene);
            d_scaled = (d+1)/2;
            ms_vis = 10 + 40*(d_scaled.^2);
            text(0.02,0.92,'\bullet','Units','normalized','Color',[0.5 0.5 0.5],'FontSize',ms_vis);
            
            % Axes formatting
            title(var_names{idx_gene},'Interpreter','tex');
            set(gca,'LineWidth',1.5,'FontSize',14,'XTick',[0 24 48],'TickDir','both',...
                'TickLength',[0.01 0.025],'Box','off'); xlim([-1,50]);
            if gi <= (nRows-1)*nCols, set(gca,'XTickLabel',[]); end
        end
        
        % Add marker size legend 
        annotation('textbox',[0.8 0.1 0.15 0.15], ...
                   'String','Marker size = correlation with cluster mean', ...
                   'EdgeColor','none','FontSize',12);

        % Add example bullets for correlations -0.5, 0, 0.5, 1.0
        ref_vals = [-0.5, 0, 0.5, 1.0];
        for i = 1:numel(ref_vals)
            d_scaled = (ref_vals(i)+1)/2;
            ms_vis = 10 + 40*(d_scaled.^2);   % same mapping as data
            text(0.82,0.12+0.06*i,'\bullet','Units','normalized', ...
                 'Color',[0.5 0.5 0.5],'FontSize',ms_vis);
            text(0.86,0.12+0.06*i,sprintf('%.1f',ref_vals(i)), ...
                 'Units','normalized','FontSize',12);
        end

        
        sgtitle(sprintf('Cluster %d',c));
    end
end

%% ==== Helper function for heatmap + dendrogram ====
function h = plotScenario(titleStr, Z, leafOrder, var_names, Zmat_timeXgene, time)
    nums = 14;
    map = zeros(2*nums, 3);
    for ind = 1:nums
        t = (ind - 1) / (nums - 1);
        map(ind, :) = [1, t, t].^3;
    end
    for ind = 1:nums
        t = (ind - 1) / (nums - 1);
        map(ind+nums,:) = [1-t,1-t,1].^3;
    end
    ordered_var_names = var_names(leafOrder);
    Zmat_ord = Zmat_timeXgene(:, leafOrder);
    figure('Name',titleStr,'Color','w','Position',[100 100 1100 800]);
    subplot(2,1,1);
    [h,~,~] = dendrogram(Z,0,'Reorder',leafOrder,'Orientation','top');
    set(h,'Color',[0 0 0],'LineWidth',2);
    title([titleStr ' - Dendrogram']);
    subplot(2,1,2);
    imagesc(1:numel(ordered_var_names), time, Zmat_ord);
    colormap(map); colorbar('LineWidth',1);
    set(gca,'CLim',[-4.5 4.5],'FontSize',14,'YDir','normal','XTick',1:numel(ordered_var_names),...
        'XTickLabel',ordered_var_names,'XTickLabelRotation',90,'TickLabelInterpreter','tex',...
        'LineWidth',1,'YTick',[0,12,24,36,48],'YTickLabel',[0,12,24,36,48]);
    ylabel('Time (h)'); 
end
