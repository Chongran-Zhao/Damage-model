function out = cross_otimes_2d_to_6d(M1, M2, M3)
out = zeros(3,3,3,3,3,3);
for ii = 1:3
    for jj = 1:3
        for kk = 1:3
            for ll = 1:3
                for mm = 1:3
                    for nn = 1:3
                        out(ii,jj,kk,ll,mm,nn) = M1(ii,jj) * M2(kk,ll) * M3(mm,nn);
                    end
                end
            end
        end
    end
end
end