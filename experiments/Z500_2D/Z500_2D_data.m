%% Daily Z500 anomalies over 30-70N (9 x 36 grid points) from the NCEP/NCAR reanalysis 1979-2024, downloaded from NOAA PSL and cached per year
function data_file = Z500_2D_data()
    d = fullfile(Repo_root, 'experiments', 'Z500_2D');
    data_file = fullfile(d, 'data', 'Z500_2D_data.mat');
    if isfile(data_file)
        return
    end
    years = 1979:2024;
    yr_train = 1979:2014;
    url_fmt = 'https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis.dailyavgs/pressure/hgt.%d.nc';
    lat_rng = 9:25;
    lev_idx = 6;
    lon_sel = 1:4:144;
    lat_sel = 1:2:17;
    cache_dir = fullfile(d, 'data', 'ncep_cache');
    if ~isfolder(cache_dir)
        mkdir(cache_dir);
    end

    Z_y = cell(numel(years), 1);
    dn_y = cell(numel(years), 1);
    for iy = 1:numel(years)
        y = years(iy);
        cache_file = fullfile(cache_dir, sprintf('z500band_%d.mat', y));
        if isfile(cache_file)
            load(cache_file, 'Zy', 'dny');
        else
            url = sprintf(url_fmt, y);
            ok = false;
            for attempt = 1:5
                try
                    hgt = ncread(url, 'hgt', [1, lat_rng(1), lev_idx, 1], [Inf, numel(lat_rng), 1, Inf]);
                    tim = ncread(url, 'time');
                    ok = true;
                    break
                catch
                    fprintf('  retry %d for %d\n', attempt, y);
                    pause(10);
                end
            end
            if ~ok
                error('download failed for %d', y);
            end
            Zy = permute(double(squeeze(hgt)), [3, 1, 2]);
            dny = datenum(1800, 1, 1) + double(tim) / 24;
            save(cache_file, 'Zy', 'dny');
        end
        Z_y{iy} = Zy;
        dn_y{iy} = dny(:);
        fprintf('year %d: %d days\n', y, size(Zy, 1));
    end
    Z_all = cat(1, Z_y{:});
    dn_all = vertcat(dn_y{:});

    dv = datevec(dn_all);
    keep = ~(dv(:, 2) == 2 & dv(:, 3) == 29);
    Z = Z_all(keep, lon_sel, lat_sel);
    dn = dn_all(keep);
    lon = (lon_sel - 1) * 2.5;
    lat = 90 - (lat_rng(lat_sel) - 1) * 2.5;
    Nlon = numel(lon);
    Nlat = numel(lat);
    Z = reshape(Z, size(Z, 1), Nlon * Nlat);
    dv = datevec(dn);
    lp = (mod(dv(:, 1), 4) == 0 & mod(dv(:, 1), 100) ~= 0) | mod(dv(:, 1), 400) == 0;
    doy = floor(dn) - datenum(dv(:, 1), 1, 1) + 1;
    doy(lp & dv(:, 2) > 2) = doy(lp & dv(:, 2) > 2) - 1;
    in_train = ismember(dv(:, 1), yr_train);
    A = zeros(size(Z));
    for q = 1:365
        sel = doy == q;
        ref = sel & in_train;
        A(sel, :) = (Z(sel, :) - mean(Z(ref, :), 1)) ./ std(Z(ref, :), 0, 1);
    end
    T = size(A, 1);
    tt = (dn - dn(1)) / 365.25;
    G = [ones(T, 1), tt];
    coef = G(in_train, :) \ A(in_train, :);
    S = A - G * coef;
    n_tr = sum(in_train);
    n_fit = round(0.8 * n_tr);
    i_val = (n_fit + 1):n_tr;
    i_test = (n_tr + 1):T;
    save(data_file, 'S', 'lon', 'lat', 'Nlon', 'Nlat', 'dn', 'n_fit', 'i_val', 'i_test');
    fprintf('Z500 2D: %d x %d = %d grid points, %d fit / %d validation / %d test days\n', Nlat, Nlon, Nlon * Nlat, n_fit, numel(i_val), numel(i_test));
end
