function out = get_zeta_t(zeta_infty, iota, alpha_t)
out = zeros(length(alpha_t), 1);
for ii = 1:length(out)
    out(ii) = zeta_infty * (1.0 - exp(-alpha_t(ii)/iota));
end
end