function out = cross_otimes_1d_to_6d(N1, N2, N3, N4, N5, N6)
out = zeros(3,3,3,3,3,3);
for ii = 1:3
    for jj = 1:3
        for kk = 1:3
            for ll = 1:3
                for mm = 1:3
                    for nn = 1:3
                        out(ii,jj,kk,ll,mm,nn) = N1(ii) * N2(jj) * N3(kk) * N4(ll) * N5(mm) * N6(nn);
                    end
                end
            end
        end
    end
end
end