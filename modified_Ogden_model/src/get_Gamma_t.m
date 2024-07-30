function out = get_Gamma_t(mu_neq, m_neq, n_neq, eta_d, Ft, time)
out = zeros(3,3,length(time));
for ii = 1:length(time)
    tol_r = 1.0;
    tol_a = 1.0;
    counter = 1;
    if ii == 1
        C_old = eye(3);
        F_new = Ft(:,:,ii);
        C_new = F_new' * F_new;
        C_mid = 0.5 .* (C_old + C_new);
        Gamma_old = eye(3);
        Gamma_new = eye(3);
        dt = time(ii);
    else
        F_old = Ft(:,:,ii-1);
        F_new = Ft(:,:,ii);
        C_old = F_old' * F_old;
        C_new = F_new' * F_new;
        C_mid = 0.5 .* (C_old + C_new);
        Gamma_old = out(:,:,ii-1);
        Gamma_new = out(:,:,ii-1);
        dt = time(ii) - time(ii-1);
    end
    Gamma_mid = 0.5 .* (Gamma_old + Gamma_new);
    tangent0 = get_res_tangent(mu_neq, m_neq, n_neq, eta_d, C_mid, Gamma_mid, dt);
    residual0 = get_res(mu_neq, m_neq, n_neq, eta_d, C_mid, Gamma_old, Gamma_new, dt);
    incremental0 = solve_AB(tangent0, -residual0);
    Gamma_new = Gamma_new + incremental0;
    while (tol_a > 1e-10 && counter < 10) || (tol_r > 1e-10 && counter <10)
        Gamma_mid = 0.5 .* (Gamma_old + Gamma_new);
        tangent = get_res_tangent(mu_neq, m_neq, n_neq, eta_d, C_mid, Gamma_mid, dt);
        residual = get_res(mu_neq, m_neq, n_neq, eta_d, C_mid, Gamma_old, Gamma_new, dt);
        incremental = solve_AB(tangent, -residual);
        Gamma_new = Gamma_new + incremental;
        tol_r = norm(residual) / norm(residual0);
        tol_a = norm(residual);
        counter = counter + 1;
    end
    out(:,:,ii) = Gamma_new;
    % format long
        % Gamma_old
        % disp(Gamma_new);
    % [V, D] = eig(Gamma_new)
    % disp(eta_d);
    % disp(ii)
    % if ii == 1
    %     dt = time(ii);
    %     Gamma_old = eye(3);
    %     C_old = eye(3);
    %     Q_proj = get_proj_Q(m_neq, n_neq, Gamma_old);
    %     T_neq = get_T_neq(mu_neq, m_neq, n_neq, C_old, Gamma_old);
    %     Q = contract(T_neq, Q_proj);
    %     out(:,:,ii) = Gamma_old + dt ./ eta_d .* Q;
    %     % out(:,:,ii) = eye(3);
    % else
    %     dt = time(ii) - time(ii-1);
    %     Gamma_old = out(:,:,ii-1);
    %     C_old = Ft(:,:,ii-1)' * Ft(:,:,ii-1);
    %     Q_proj = get_proj_Q(m_neq, n_neq, Gamma_old);
    %     T_neq = get_T_neq(mu_neq, m_neq, n_neq, C_old, Gamma_old);
    %     Q = contract(T_neq, Q_proj);
    %     out(:,:,ii) = Gamma_old + dt ./ eta_d .* Q;
    % end
end
end