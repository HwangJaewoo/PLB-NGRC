%% Local states (3 x 3 grid points over k days) of all grid points at the times tt
function L = Z500_2D_local_states(Saug, tt, NB, k)
    nt = numel(tt);
    [P, nnb] = size(NB);
    L = zeros(nt * P, nnb * k);
    c = 0;
    for q = 0:k - 1
        blk = Saug(tt - q, :);
        for n = 1:nnb
            c = c + 1;
            L(:, c) = reshape(blk(:, NB(:, n)), [], 1);
        end
    end
end
