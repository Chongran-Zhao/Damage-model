function out = get_S_iso_neq(mu_neq, m_neq, n_neq, C, Ev)
out = contract( get_T_neq(mu_neq, m_neq, n_neq, C, Ev), get_proj_Q(m_neq, n_neq, C));
end