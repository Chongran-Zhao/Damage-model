% solve AX = B,
% where A is the forth-order tensor with minor symmetry
% X and B should be symmetric matrix
function out = solve_AB(A, B)
% A = real(A);
% B = real(B);
map = [1, 4, 5, 4, 2, 6, 5, 6, 3];
AA = zeros(6,6);
BB = zeros(6,1);
out = zeros(3,3);
for ii = 1:3
    for jj = 1:3
        BB(map(3*ii+jj-3), 1) = B(ii,jj);
        for kk = 1:3
            for ll = 1:3
                AA(map(3*ii+jj-3), map(3*kk+ll-3)) = A(ii,jj,kk,ll);
            end
        end
    end
end

XX = AA \ BB;

for ii = 1:3
    for jj = 1:3
        if ii == jj
            out(ii,jj) = XX(map(3*ii+jj-3));
        else
            out(ii,jj) = 0.5 * XX(map(3*ii+jj-3));
        end
    end
end
end