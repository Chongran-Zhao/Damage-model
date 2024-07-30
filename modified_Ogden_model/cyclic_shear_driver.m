clc; clear; close all

addpath("src");
data = readmatrix('../exp_data_cyclic_shear_raw/400.xlsx');
time = data(:,1);
P_exp = data(:,3);
gamma = data(:,4);

% Uniaxial tension experimental case
Ft = zeros(3,3,length(time));
Ft(1,1,:) = 1.0;
Ft(2,2,:) = 1.0;
Ft(3,3,:) = 1.0;
Ft(1,2,:) = gamma(:);

% parameters
mu_eq = [1.0, 1.0];
m_eq = [1.0, 1.0];
n_eq = [1.0, 1.0];

mu_neq = [1.0, 1.0];
m_neq = [1.0, 1.0];
n_neq = [1.0, 1.0];
eta_d = [1.0, 1.0];

m = [1.0];
r = [100.0];
beta = [1.0];
[paras0, num_eq, num_neq, lb, ub] = array_to_paras(mu_eq, m_eq, n_eq, mu_neq, m_neq, n_neq, eta_d, m, r, beta);
objectiveFunction = @(paras) objective(paras, Ft, P_exp, time, num_eq, num_neq);
options = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'MaxIterations', 1000, ...
    'TolFun', 1e-10, ...
    'TolX', 1e-10, ...
    'Display', 'iter-detailed', ...
    'MaxFunctionEvaluations', 3500,  ...
    'PlotFcn', 'optimplotfval');
[paras, ~] = lsqnonlin( objectiveFunction, paras0, lb, ub, options);
plot_result(paras, num_eq, num_neq, Ft, time, P_exp);
print(gcf, '-djpeg', 'fig_cyclic_shear.jpg');