function out = get_P_iso_neq_list(mu_neq, m_neq, n_neq, eta_d, Ft, time)
out = zeros(size(Ft));
for ii = 1:length(mu_neq)
    Ev_t = get_Ev_t(mu_neq(ii), m_neq(ii), n_neq(ii), eta_d(ii), Ft, time);
    for jj = 1:length(time)
        out(:,:,jj) = out(:,:,jj) + get_P_iso_neq(mu_neq(ii), m_neq(ii), n_neq(ii), Ft(:,:,jj), Ev_t(:,:,jj));
    end
end
end
% EOF