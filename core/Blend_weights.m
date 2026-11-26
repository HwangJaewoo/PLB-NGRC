%% Normalized Gaussian blending weights of the rows of U for the balls with centers C and widths sig
function Pi = Blend_weights(U, C, sig)
    n = size(U, 1);
    B = size(C, 1);
    E = zeros(n, B);
    for b = 1:B
        E(:, b) = sum((U - C(b, :)).^2, 2) / (2 * sig(b)^2);
    end
    E = E - min(E, [], 2);
    Pi = exp(-E);
    Pi = Pi ./ sum(Pi, 2);
end
