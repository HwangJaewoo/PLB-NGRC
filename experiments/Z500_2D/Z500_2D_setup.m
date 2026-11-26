%% Data, neighbors (3 x 3 grid points), scaled latitude and seasonal phase for the Z500 data over the latitude band
function cfg = Z500_2D_setup()
    load(Z500_2D_data(), 'S', 'Nlon', 'Nlat', 'dn', 'n_fit', 'i_val', 'i_test');
    P = Nlon * Nlat;
    ell = 1;
    NB = Neighbors(Nlon, Nlat, ell);
    latP = (ceil((1:P)' / Nlon) - (Nlat + 1) / 2) / ((Nlat - 1) / 2);
    dv = datevec(dn);
    lp = (mod(dv(:, 1), 4) == 0 & mod(dv(:, 1), 100) ~= 0) | mod(dv(:, 1), 400) == 0;
    doy = floor(dn) - datenum(dv(:, 1), 1, 1) + 1;
    doy(lp & dv(:, 2) > 2) = doy(lp & dv(:, 2) > 2) - 1;
    ph = [cos(2 * pi * (doy - 1) / 365), sin(2 * pi * (doy - 1) / 365)];
    cfg = struct('S', S, 'Saug', [S, zeros(size(S, 1), 1)], 'Nlon', Nlon, 'Nlat', Nlat, 'P', P, ...
        'NB', NB, 'latP', latP, 'ph', ph, 'n_fit', n_fit, 'i_val', i_val, 'i_test', i_test, ...
        'div_tol', 50, 'chunk', 200, 'n_sub', 2e5, 'max_lag', 10);
end

function NB = Neighbors(Nlon, Nlat, ell)
    P = Nlon * Nlat;
    NB = zeros(P, (2 * ell + 1)^2);
    for ilat = 1:Nlat
        for ilon = 1:Nlon
            p = ilon + Nlon * (ilat - 1);
            n = 0;
            for a = -ell:ell
                for b = -ell:ell
                    n = n + 1;
                    jl = ilat + a;
                    jo = mod(ilon + b - 1, Nlon) + 1;
                    if jl < 1 || jl > Nlat
                        NB(p, n) = P + 1;
                    else
                        NB(p, n) = jo + Nlon * (jl - 1);
                    end
                end
            end
        end
    end
end
