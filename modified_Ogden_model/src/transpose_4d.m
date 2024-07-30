function out = transpose_4d(A)
out = zeros(3,3,3,3);
for ii = 1:3
    for jj = 1:3
        for kk = 1:3
            for ll = 1:3
                out(ii,jj,kk,ll) = A(kk,ll,ii,jj);
            end
        end
    end
end
end