function out = get_eta_t(mu_eq, m_eq, n_eq, mu_neq, m_neq, n_neq, eta_d, Ft, time, m, r, beta)
out = ones(length(time), 1);
Psi_hat = zeros(length(time), 1);
max_Psi_hat = 0.0;
for ii = 1:length(mu_eq)
    for jj = 1:length(time)
        F = Ft(:,:,jj);
        C = F' * F;
        Psi_hat(jj) = Psi_hat(jj) + get_CR_energy_eq(mu_eq(ii), m_eq(ii), n_eq(ii), C);
    end
end

for ii = 1:length(mu_neq)
    Ev_t = get_Ev_t(mu_neq(ii), m_neq(ii), n_neq(ii), eta_d(ii), Ft, time);
    for jj = 1:length(time)
        F = Ft(:,:,jj);
        C = F' * F;
        Psi_hat(jj) = Psi_hat(jj) + get_CR_energy_neq(mu_neq(ii), m_neq(ii), n_neq(ii), C, Ev_t(:,:,jj));
    end
end
for ii = 1:length(time)
    if Psi_hat(ii) >= max_Psi_hat
        max_Psi_hat = Psi_hat(ii);
    end
    out(ii) = 1.0 - erf( (max_Psi_hat -  Psi_hat(ii)) / (m +  beta * max_Psi_hat) ) / r;
end
end