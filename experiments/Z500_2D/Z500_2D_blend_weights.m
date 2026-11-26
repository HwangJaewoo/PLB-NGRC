%% Normalized Gaussian blending weights of the conditioned local states Z for the balls with centers C and widths sig
function w = Z500_2D_blend_weights(L, C, sig)
    if size(C, 1) == 1
        w = ones(size(L, 1), 1);
        return
    end
    lg = -pdist2(L, C).^2 ./ (2 * sig(:)'.^2);
    lg = lg - max(lg, [], 2);
    g = exp(lg);
    w = g ./ sum(g, 2);
end
