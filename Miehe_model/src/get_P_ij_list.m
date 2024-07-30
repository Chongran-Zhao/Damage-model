function out = get_P_ij_list(ii, jj, mu_eq, m_eq, n_eq, mu_neq, m_neq, n_neq, eta_d, d_a, eta_a, d_b, eta_b, Ft, time)
out = zeros(length(time), 1);
P_pre = get_P_list(mu_eq, m_eq, n_eq, mu_neq, m_neq, n_neq, eta_d, Ft, time);
d_t = get_d_t(mu_eq, m_eq, n_eq, mu_neq, m_neq, n_neq, eta_d, Ft, time, d_a, eta_a, d_b, eta_b);
for kk = 1:length(time)
out(kk) = (1.0 - d_t(kk)) * P_pre(ii,jj,kk);
end
end