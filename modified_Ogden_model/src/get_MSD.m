function out = get_MSD(P1_exp, P1_list)
out = 0.0;
for ii = 1:length(P1_exp)
    out = out + (P1_exp(ii)-P1_list(ii)) .* (P1_exp(ii)-P1_list(ii))/ length(P1_exp);
end
end