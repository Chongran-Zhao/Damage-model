% refer to eqn. 2.8 of Liu, Guan, Zhao & Luo 2024 preprint
% The input C can also be replaced by Gamma, because the transformation
% of E and C is same with Ev and Gamma
function out = get_proj_Q(m, n, C)
out = zeros(3,3,3,3);
[V, D] = eig(C);
N1 = V(:,1);
N2 = V(:,2);
N3 = V(:,3);
lambda1 = sqrt(D(1,1));
lambda2 = sqrt(D(2,2));
lambda3 = sqrt(D(3,3));

E1 = get_CR_scale(m, n, lambda1);
E2 = get_CR_scale(m, n, lambda2);
E3 = get_CR_scale(m, n, lambda3);

E1_der_1 = get_CR_scale_der(m, n, lambda1, 1);
E2_der_1 = get_CR_scale_der(m, n, lambda2, 1);
E3_der_1 = get_CR_scale_der(m, n, lambda3, 1);

d1 = E1_der_1 / lambda1;
d2 = E2_der_1 / lambda2;
d3 = E3_der_1 / lambda3;

if lambda1 == lambda2
    theta_12 = d1;
else
    theta_12 = 2.0 * (E1 - E2) / (lambda1*lambda1 - lambda2*lambda2);
end

if lambda2 == lambda3
    theta_23 = d2;
else
    theta_23 = 2.0 * (E2 - E3) / (lambda2*lambda2 - lambda3*lambda3);
end

if lambda1 == lambda3
    theta_13 = d1;
else
    theta_13 = 2.0 * (E1 - E3) / (lambda1*lambda1 - lambda3*lambda3);
end

out = out + d1 .* cross_otimes_1d_to_4d(N1, N1, N1, N1);
out = out + d2 .* cross_otimes_1d_to_4d(N2, N2, N2, N2);
out = out + d3 .* cross_otimes_1d_to_4d(N3, N3, N3, N3);

out = out + theta_12 .* dot_otimes(N1, N2) + theta_12 .* dot_otimes(N2, N1);
out = out + theta_13 .* dot_otimes(N1, N3) + theta_13 .* dot_otimes(N3, N1);
out = out + theta_23 .* dot_otimes(N2, N3) + theta_23 .* dot_otimes(N3, N2);
% if m == 0.0 && n == 0.0
%     out = 0.5 .* (get_proj_sub_Q(0, C) + get_proj_sub_Q(0, C));
% elseif m == 0 && n ~= 0
%     out = get_proj_sub_Q(-n, C);
% elseif m ~= 0 && n == 0
%     out = get_proj_sub_Q(m, C);
% else
%     out = (m .* get_proj_sub_Q(m, C) ./ (m+n)) + (n .* get_proj_sub_Q(-n, C) ./ (m+n));
% end
end