function out = get_d_t(mu_eq, m_eq, n_eq, mu_neq, m_neq, n_neq, eta_d, Ft, time, d_a, eta_a, d_b, eta_b)
out = ones(length(time), 1);
Psi = zeros(length(time), 1);
beta_t = zeros(length(time), 1);
for ii = 1:length(mu_eq)
    for jj = 1:length(time)
        F = Ft(:,:,jj);
        C = F' * F;
        Psi(jj) = Psi(jj) + get_CR_energy_eq(mu_eq(ii), m_eq(ii), n_eq(ii), C);
    end
end

for ii = 1:length(mu_neq)
    Ev_t = get_Ev_t(mu_neq(ii), m_neq(ii), n_neq(ii), eta_d(ii), Ft, time);
    for jj = 1:length(time)
        F = Ft(:,:,jj);
        C = F' * F;
        Psi(jj) = Psi(jj) + get_CR_energy_neq(mu_neq(ii), m_neq(ii), n_neq(ii), C, Ev_t(:,:,jj));
    end
end
max_Psi = Psi(1);
beta_t(1) = Psi(1);
for ii = 2:length(time)
    if Psi(ii) >= max_Psi
        max_Psi = Psi(ii);
    end
    beta_t(ii) = Psi(ii-1) + abs(Psi(ii)-Psi(ii-1));
end
for ii = 1:length(time)
    out(ii) = d_a * (1.0 - exp(-max_Psi / eta_a)) + d_b * (1.0 - exp(-beta_t(ii) / eta_b));
end
end