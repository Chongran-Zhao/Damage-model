function out = get_P_iso_eq_list(mu_eq_list, m_eq_list, n_eq_list, Ft)
out = zeros(size(Ft));
if length(mu_eq_list) == length(m_eq_list) && length(m_eq_list) == length(n_eq_list)
    for ii = 1:length(mu_eq_list)
        for jj = 1:size(Ft, 3)
            out(:,:,jj) = out(:,:,jj) + get_P_iso_eq(mu_eq_list(ii), m_eq_list(ii), n_eq_list(ii), Ft(:,:,jj));
        end
    end
else
    error("ERROR: Length of input parameters of equilibrium part don't match!");
end

end
% EOF