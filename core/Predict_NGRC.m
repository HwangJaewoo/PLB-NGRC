%% One-step prediction of NGRC (PLB-NGRC, P-NGRC for s = 0, LB-NGRC on the state of all sites) for the local states U
function Yhat = Predict_NGRC(mdl, U)
    ord = 2;
    if isfield(mdl, 'ord')
        ord = mdl.ord;
    end
    Pi = Blend_weights(U, mdl.C, mdl.sig);
    Yhat = zeros(size(U, 1), mdl.m);
    for b = 1:size(mdl.C, 1)
        Yhat = Yhat + Pi(:, b) .* (Poly_features(U - mdl.C(b, :), ord) * mdl.W(:, :, b));
    end
end
