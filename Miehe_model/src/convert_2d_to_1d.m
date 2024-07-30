function out = convert_2d_to_1d(A)
out = zeros(9,1);
for ii = 1:3
    for jj = 1:3
        out(3*(ii-1)+jj) = A(ii,jj);
    end
end
end