% The function is used to obtain six-order projection tensor L
% refer to eqn. 2.11 of Liu, Guan, Zhao & Luo 2024 preprint
function out = get_proj_L(m, n, C)
out = zeros(3,3,3,3,3,3);

[V, D] = eig(C);
lambda1 = sqrt(D(1,1));
lambda2 = sqrt(D(2,2));
lambda3 = sqrt(D(3,3));

N1 = V(:,1);
N2 = V(:,2);
N3 = V(:,3);

E1 = get_CR_scale(m, n, lambda1);
E2 = get_CR_scale(m, n, lambda2);
E3 = get_CR_scale(m, n, lambda3);

E1_der_1 = get_CR_scale_der(m, n, lambda1, 1);
E2_der_1 = get_CR_scale_der(m, n, lambda2, 1);
E3_der_1 = get_CR_scale_der(m, n, lambda3, 1);

E1_der_2 = get_CR_scale_der(m, n, lambda1, 2);
E2_der_2 = get_CR_scale_der(m, n, lambda2, 2);
E3_der_2 = get_CR_scale_der(m, n, lambda3, 2);

% first part
f1 = -E1_der_1 / (lambda1^3) + E1_der_2 / (lambda1^2);
f2 = -E2_der_1 / (lambda2^3) + E2_der_2 / (lambda2^2);
f3 = -E3_der_1 / (lambda3^3) + E3_der_2 / (lambda3^2);

out = out + f1 .* cross_otimes_1d_to_6d(N1, N1, N1, N1, N1, N1);
out = out + f2 .* cross_otimes_1d_to_6d(N2, N2, N2, N2, N2, N2);
out = out + f3 .* cross_otimes_1d_to_6d(N3, N3, N3, N3, N3, N3);

% get d
d1 = E1_der_1 / lambda1;
d2 = E2_der_1 / lambda2;
d3 = E3_der_1 / lambda3;

% get theta and xi
if lambda1 == lambda2
    theta_12 = d1;
    theta_21 = d2;
    xi_12 = 0.125 * f1;
    xi_21 = 0.125 * f2;
else
    theta_12 = 2.0 * (E1 - E2) / (lambda1*lambda1 - lambda2*lambda2);
    theta_21 = theta_12;
    xi_12 = 0.5 * (theta_12 - d2) / (lambda1*lambda1 - lambda2*lambda2);
    xi_21 = 0.5 * (theta_21 - d1) / (lambda2*lambda2 - lambda1*lambda1);
end

if lambda2 == lambda3
    theta_23 = d2;
    theta_32 = d3;
    xi_23 = 0.125 * f2;
    xi_32 = 0.125 * f3;
else
    theta_23 = 2.0 * (E2 - E3) / (lambda2*lambda2 - lambda3*lambda3);
    theta_32 = theta_23;
    xi_23 = 0.5 * (theta_23 - d3) / (lambda2*lambda2 - lambda3*lambda3);
    xi_32 = 0.5 * (theta_32 - d2) / (lambda3*lambda3 - lambda2*lambda2);
end

if lambda1 == lambda3
    theta_13 = d1;
    theta_31 = d3;
    xi_13 = 0.125 * d1;
    xi_31 = 0.125 * d3;
else
    theta_13 = 2.0 * (E1 - E3) / (lambda1*lambda1 - lambda3*lambda3);
    theta_31 = theta_13;
    xi_13 = 0.5 * (theta_13 - d3) / (lambda1*lambda1 - lambda3*lambda3);
    xi_31 = 0.5 * (theta_31 - d1) / (lambda3*lambda3 - lambda1*lambda1);
end

% second part
out = out + xi_12 .* (get_H_abc(N1, N2, N2) + get_H_abc(N2, N1, N2) + get_H_abc(N2, N2, N1));
out = out + xi_13 .* (get_H_abc(N1, N3, N3) + get_H_abc(N3, N1, N3) + get_H_abc(N3, N3, N1));

out = out + xi_21 .* (get_H_abc(N2, N1, N1) + get_H_abc(N1, N2, N1) + get_H_abc(N1, N1, N2));
out = out + xi_23 .* (get_H_abc(N2, N3, N3) + get_H_abc(N3, N2, N3) + get_H_abc(N3, N3, N2));

out = out + xi_31 .* (get_H_abc(N3, N1, N1) + get_H_abc(N1, N3, N1) + get_H_abc(N1, N1, N3));
out = out + xi_32 .* (get_H_abc(N3, N2, N2) + get_H_abc(N2, N3, N2) + get_H_abc(N2, N2, N3));

% get eta
if lambda1 == lambda2 && lambda2 == lambda3
    eta = 0.125 * f1;
end

if lambda1 == lambda2
    eta = xi_31;
elseif lambda2 == lambda3
    eta = xi_12;
elseif lambda3 == lambda1
    eta = xi_23;
else
    eta = 0.5 * E1 / ((lambda1*lambda1 - lambda2*lambda2) * (lambda1*lambda1 - lambda3*lambda3))...
        + 0.5 * E2 / ((lambda2*lambda2 - lambda1*lambda1) * (lambda2*lambda2 - lambda3*lambda3))...
        + 0.5 * E3 / ((lambda3*lambda3 - lambda1*lambda1) * (lambda3*lambda3 - lambda2*lambda2));
end

out = out + eta .* (get_H_abc(N1, N2, N3) + get_H_abc(N1, N3, N2));
out = out + eta .* (get_H_abc(N2, N1, N3) + get_H_abc(N2, N3, N1));
out = out + eta .* (get_H_abc(N3, N1, N2) + get_H_abc(N3, N2, N1));
end