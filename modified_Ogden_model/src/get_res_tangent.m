% refer to Local Newton-Raphson iteration box of Liu, Guan, Zhao & Luo 2024 preprint
function out = get_res_tangent(mu_neq, m_neq, n_neq, eta_d, C_mid, Gamma_mid, dt)
out = get_sym_idn_4d();
Q_proj = get_proj_Q(m_neq, n_neq, Gamma_mid);
T = get_T_neq(mu_neq, m_neq, n_neq, C_mid, Gamma_mid);
L_proj = get_proj_L(m_neq, n_neq, Gamma_mid);
K_proj = -2.0 .* mu_neq .* contract(transpose_4d(Q_proj), Q_proj) + contract(T, L_proj);
out = out - 0.25 .* dt / eta_d .* K_proj;
end