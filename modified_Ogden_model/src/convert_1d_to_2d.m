function out = convert_1d_to_2d(A)
out = zeros(3,3);
for ii = 1:3
    for jj = 1:3
        out(ii,jj) = A(3*(ii-1)+jj);
    end
end
end