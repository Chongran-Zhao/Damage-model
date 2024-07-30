function out = convert_4d_to_2d(A)
out = zeros(9,9);
for ii = 1:3
    for jj = 1:3
        for kk = 1:3
            for ll = 1:3
                out(3*(ii-1)+jj, 3*(kk-1)+ll) = A(ii,jj,kk,ll);
            end
        end
    end
end
end