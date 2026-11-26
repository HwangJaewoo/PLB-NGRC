%% Local states of all sites: states of the 2*nn+1 neighboring sites over k delays, one row per site and time
function U = Delay_input(S, nn, m, k)
    T = size(S, 1);
    N = size(S, 2) / m;
    d0 = (2 * nn + 1) * m;
    U = zeros((T - k + 1) * N, d0 * k);
    for j = 0:k - 1
        U(:, j * d0 + (1:d0)) = Local_input(S((k - j):(T - j), :), nn, m);
    end
end

function U = Local_input(S, nn, m)
    [T, Nm] = size(S);
    N = Nm / m;
    offs = -nn:nn;
    d0 = numel(offs) * m;
    U3 = zeros(T, N, d0);
    col = 0;
    for j = 1:numel(offs)
        for v = 1:m
            col = col + 1;
            U3(:, :, col) = circshift(S(:, (v - 1) * N + (1:N)), [0, -offs(j)]);
        end
    end
    U = reshape(U3, T * N, d0);
end
