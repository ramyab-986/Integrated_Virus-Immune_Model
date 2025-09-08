clear; close all; clc;

load('param.mat');
load('SteadyState_120h.mat', 'Tss', 'Yss');

y0 = Yss(end,:);
y0(2) = 100;           % Virus Input  

tend = 200*60;
tspan = linspace(0,tend);

% Parameter ranges (log-spaced)
c = linspace(log10(0.03), log10(30), 31);   
VC_values = 10.^c;                          % Viral antagonism values
Ia_values = 10.^c;                          % Immune activation values
I_n = 1;

variable_values = zeros(length(Ia_values), length(VC_values));

parfor i = 1:length(Ia_values)  % Outer loop over I_a
    I_a = Ia_values(i);
    
    local_results = zeros(1, length(VC_values)); 
    
    for j = 1:length(VC_values) % Inner loop over VC
        VC = VC_values(j);
        
        % Run ODE solver with current parameters
        [~, Y] = ode23s(@(t,y) ODEs(t, y, param, I_n, I_a, VC), tspan, y0);
        
        % Store the Viral load values
        local_results(j) = Y(end, 1);
    end
    
    variable_values(i,:) = local_results; 
end
save('variable_values_31.mat', 'VC_values', 'Ia_values', 'variable_values');

%% Plotting

load('variable_values_31.mat', 'VC_values', 'Ia_values', 'variable_values');

% Create a heatmap plot
figure;
imagesc(log10(Ia_values), log10(VC_values), log10(variable_values));
colorbar;
set(gca, 'FontSize', 24, 'LineWidth', 2)
xlabel('log_{10}(V_{I, N} values)');
ylabel('log_{10}(I_{V, N} values)');
title('HCV')

% Set axis properties
set(gca, 'YDir', 'normal');  
colormap('sky');  
c = colorbar;
clim([-2 3]);  
c.Ticks = [-2, -1, 0, 1, 2, 3];
c.TickLabels = arrayfun(@(x) ['10^{' num2str(x) '}'], c.Ticks, 'UniformOutput', false);


