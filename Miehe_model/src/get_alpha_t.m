function out = get_alpha_t(mu_eq, m_eq, n_eq, mu_neq, m_neq, n_neq, eta_d, Ft, time)
out = zeros(length(time), 1);
Psi = zeros(length(time), 1);
max_value = 0.0;
for ii = 1:length(time)
    F = Ft(:,:,ii);
    C = F' * F;
    for jj = 1:length(mu_eq)
        Psi(ii) = Psi(ii) + get_CR_energy_eq(mu_eq(jj), m_eq(jj), n_eq(jj), C);
    end
    for jj = 1:length(mu_neq)
        Ev_t = get_Ev_t(mu_neq, m_neq, n_neq, eta_d(jj), Ft, time);
        Psi(ii) = Psi(ii) + get_CR_energy_neq(mu_neq(jj), m_neq(jj), n_neq(jj), C, Ev_t(:,:,ii));
    end
    if Psi(ii) >= max_value
        max_value = Psi(ii);
        out(ii) = max_value;
    else
        out(ii) = max_value;
    end
end
end