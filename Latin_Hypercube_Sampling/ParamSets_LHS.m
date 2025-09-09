clear; clc; close all; 

iter_prm = 25000; 
iter_ini = 1;
load('param.mat', 'param');

lb_prm = param;
ub_prm = param;

arr_ch = [3:110, 112:129];
lb_prm(1, arr_ch) = 0.1 * param(1, arr_ch);
ub_prm(1, arr_ch) = 10 * param(1, arr_ch);

y0 = zeros(1,75);       
y0(2) = 100;           % virus input
y0(10) = 5.3422;
y0(12) = 277.78;
y0(14) = 3.0832;
y0(16) = 97.169;
y0(18) = 37.862;
y0(20) = 37.969;
y0(22) = 11.3587;
y0(25) = 101.735;
y0(27) = 24.92;
y0(34) = 151.88;
y0(36) = 1114.8;
y0(37) = 20;
y0(40) = 6.5019;
y0(41) = 20.701;
y0(44) = 10^3;
y0(45) = 10^3;
y0(47) = 45;
y0(55) = 40;
y0(58) = 41.96;
y0(63) = 500; 

lb_ini = y0*1;
ub_ini = y0*1;

X_prm = lhsdesign(iter_prm, size(lb_prm, 2),'iterations', 5);
X_ini = lhsdesign(iter_ini, size(lb_ini, 2),'iterations', 5);


for ind_prm = 1:iter_prm
    disp(ind_prm)
    for id_prm = 1:size(lb_prm, 2)
        PARAMETER(ind_prm, id_prm) = lb_prm(1, id_prm) + ...
            (ub_prm(1, id_prm) - lb_prm(1, id_prm)) * X_prm(ind_prm, id_prm);
    end
end

for ind_ini = 1:iter_ini
    for id_ini = 1:size(lb_ini, 2)
        INI_COND(ind_ini, id_ini) = lb_ini(1, id_ini) + ...
            (ub_ini(1, id_ini) - lb_ini(1, id_ini)) * X_ini(ind_ini, id_ini);
    end
end

save('bounds.mat');
