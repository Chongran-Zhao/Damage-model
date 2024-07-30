function out = get_Ev_t(mu_neq, m_neq, n_neq, eta_d, Ft, time)
out = zeros(size(Ft));
% backward Euler
% for ii = 1:length(time)
%     if ii == 1
%         dt = time(1);
%     else
%         dt = time(ii) - time(ii-1);
%     end
%     C_new = Ft(:,:,ii)' * Ft(:,:,ii);
%     E_new = get_CR_strain(m_neq, n_neq, C_new);
%     out(:,:,ii) = (mu_neq * dt / (eta_d + mu_neq * dt)) .* E_new;
%     if ii == 1
%         out(:,:,ii) = out(:,:,ii) + (eta_d / (eta_d + mu_neq * dt)) .* zeros(3,3);
%     else
%         out(:,:,ii) = out(:,:,ii) + (eta_d / (eta_d + mu_neq * dt)) .* out(:,:,ii-1);
%     end
% end

% forward Euler
% for ii = 1:length(time)
%     if ii == 1
%         dt = time(1);
%         Ev_old = zeros(3,3);
%         C_old = eye(3);
%     else
%         dt = time(ii) - time(ii-1);
%         Ev_old = out(:,:,ii-1);
%         F_old = Ft(:,:,ii-1);
%         C_old = F_old' * F_old;
%     end
%     E_old = get_CR_strain(m_neq, n_neq, C_old);
%     out(:,:,ii) = mu_neq / eta_d * dt .* (E_old - Ev_old) + Ev_old;
% end

% exponential integration
% for ii = 1:length(time)
%     if ii == 1
%         dt = time(1);
%         C_old = eye(3);
%         Ev_old = zeros(3,3);
%     else
%         dt = time(ii) - time(ii-1);
%         F_old = Ft(:,:,ii-1);
%         C_old = F_old' * F_old;
%         Ev_old = out(:,:,ii-1);
%     end
%     E_old = get_CR_strain(m_neq, n_neq, C_old);
%     out(:,:,ii) = exp(-mu_neq*dt/eta_d) .* Ev_old + exp(-mu_neq*dt/eta_d) * mu_neq / eta_d * dt .* E_old;
% end

% exponential integration
% for ii = 1:length(time)
%     if ii == 1
%         dt = time(1);
%         F_new = Ft(:,:,ii);
%         C_new = F_new' * F_new;
%         Ev_old = zeros(3,3);
%     else
%         dt = time(ii) - time(ii-1);
%         F_new = Ft(:,:,ii);
%         C_new = F_new' * F_new;
%         Ev_old = out(:,:,ii-1);
%     end
%     E_new = get_CR_strain(m_neq, n_neq, C_new);
%     out(:,:,ii) = exp(-mu_neq*dt/eta_d) .* Ev_old + mu_neq / eta_d * dt .* E_new;
% end

% exponential integration
for ii = 1:length(time)
    if ii == 1
        dt = time(1);
        C_old = eye(3);
        F_new = Ft(:,:,ii);
        C_new = F_new' * F_new;
        Ev_old = zeros(3,3);
    else
        dt = time(ii) - time(ii-1);
        F_old = Ft(:,:,ii-1);
        C_old = F_old' * F_old;
        F_new = Ft(:,:,ii);
        C_new = F_new' * F_new;
        Ev_old = out(:,:,ii-1);
    end
    E_old = get_CR_strain(m_neq, n_neq, C_old);
    E_new = get_CR_strain(m_neq, n_neq, C_new);
    out(:,:,ii) = exp(-mu_neq*dt/eta_d) .* Ev_old + 0.5 * mu_neq / eta_d * dt .* ( exp(-mu_neq*dt/eta_d) .* E_old + E_new );
end
end