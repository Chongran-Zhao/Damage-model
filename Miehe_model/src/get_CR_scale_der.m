function out = get_CR_scale_der(m, n, lambda, order)
switch order
    case 1
        out = m/(m+n) * lambda^(m-1) + n/(m+n) * lambda^(-n-1);
    case 2
        out = (m-1)*m/(m+n) * lambda^(m-2) - (n+1)*n/(m+n) * lambda^(-n-2);
end
end