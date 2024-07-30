% The function is used to obtain the stress-like tensor T of equilibrium part,
% which is conjugated with generalized strain E
% refer to sec. 2.3 of Liu, Guan, Zhao & Luo 2024 preprint
function out = get_T_eq(mu_eq, m_eq, n_eq, C)
out = 2.0 .* mu_eq .* get_CR_strain(m_eq, n_eq, C);
end