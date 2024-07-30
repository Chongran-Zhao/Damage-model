function out= convert_2d_to_4d(A)
out = zeros(3,3,3,3);
for ii = 1:3
    for jj = 1:3
        for kk = 1:3
            for ll = 1:3
                out(ii,jj,kk,ll) = A(3*(ii-1)+jj, 3*(kk-1)+ll);
            end
        end
    end
end
end