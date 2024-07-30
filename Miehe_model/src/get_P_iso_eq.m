function out = get_P_iso_eq(mu_eq, m_eq, n_eq, F)
C = F' * F;
out = F * get_S_iso_eq(mu_eq, m_eq, n_eq, C);
end