clear all; clc; close all; warning off; format shortG;

%% Housekeeping
tol = 3e-2;
load('bounds.mat', 'PARAMETER', 'INI_COND');

iter_prm = size(PARAMETER, 1);
iter_ini = size(INI_COND, 1);
t_ss = 48 * 60;
tspan = [0, t_ss];
epsilon = tol;
minpts = max(floor(iter_ini / 10), 5);

%% Preallocate memory-efficient storage
SS = single(zeros(iter_ini, 75, iter_prm));

%% Parallel computation

I_n = 1; 
I_a = 1.01388; VC = 0.439;  %Critical values
 
tic
parfor ind_prm = 1:iter_prm
    prm_arr = PARAMETER(ind_prm, :);
    prm = prm_arr';

    SS_local = zeros(iter_ini, 75, 'single');

    for ind_ini = 1:iter_ini
        y0 = INI_COND(ind_ini, :);
        [t, y] = ode23s(@(t, y) ODEs(t, y, prm, I_n, I_a, VC), tspan, y0);

        y_ss = y(end, :);  % Steady state
        SS_local(ind_ini, :) = single(y_ss);

        % SS(ind_ini, :, ind_prm) = y_ss;

   end

    SS(:, :, ind_prm) = SS_local;

    fprintf('Completed PARAMETER set %d/%d\n', ind_prm, iter_prm);

    % Optional: save partial result to avoid memory overload
    % save(sprintf('SS_part_%d.mat', ind_prm), 'SS_local');
end
toc

%% Save the full result (if memory allows)
save('session.mat', 'SS', '-v7.3');  % Use -v7.3 for large files
