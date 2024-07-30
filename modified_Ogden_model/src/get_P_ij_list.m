function out = get_P_ij_list(ii, jj, mu_eq, m_eq, n_eq, mu_neq, m_neq, n_neq, eta_d, m, r, beta, Ft, time)
out = zeros(length(time), 1);
P_pre = get_P_list(mu_eq, m_eq, n_eq, mu_neq, m_neq, n_neq, eta_d, Ft, time);
eta_t = get_eta_t(mu_eq, m_eq, n_eq, mu_neq, m_neq, n_neq, eta_d, Ft, time, m, r, beta);
for kk = 1:length(time)
out(kk) = eta_t(kk) * P_pre(ii,jj,kk);
end
end