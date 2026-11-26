%% Polynomial features [1, z, quadratic monomials, (cubic monomials if ord >= 3)] of each row of Z
function O = Poly_features(Z, ord)
    [n, d] = size(Z);
    nq = d * (d + 1) / 2;
    Q = zeros(n, nq);
    c = 0;
    for i = 1:d
        for j = i:d
            c = c + 1;
            Q(:, c) = Z(:, i) .* Z(:, j);
        end
    end
    O = [ones(n, 1), Z, Q];
    if ord >= 3
        nc = d * (d + 1) * (d + 2) / 6;
        Q3 = zeros(n, nc);
        c = 0;
        for i = 1:d
            for j = i:d
                for l = j:d
                    c = c + 1;
                    Q3(:, c) = Z(:, i) .* Z(:, j) .* Z(:, l);
                end
            end
        end
        O = [O, Q3];
    end
end
