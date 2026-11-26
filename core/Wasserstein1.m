%% Wasserstein distance W1 between the distributions of A and B, normalized by the 1-99% quantile range of A, averaged over the m variables
function d = Wasserstein1(A, B, nq, m)
    if nargin < 4
        m = 1;
    end
    N = size(A, 2) / m;
    q = ((1:nq) - 0.5) / nq;
    d = 0;
    for v = 1:m
        cols = (v - 1) * N + (1:N);
        a = A(:, cols);
        b = B(:, cols);
        scale = max(diff(quantile(a(:), [0.01, 0.99])), 1e-12);
        d = d + mean(abs(quantile(a(:), q) - quantile(b(:), q))) / scale;
    end
    d = d / m;
end
