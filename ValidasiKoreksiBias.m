%% ============================================================
%  VALIDASI OUT-OF-SAMPLE KOREKSI BIAS SUHU CMIP6 PULAU JAWA
%  Periode validasi 2000-2014, fungsi transfer dari kalibrasi 1970-1999
%
%  Prinsip. Fungsi transfer QM dan QDM dilatih HANYA pada data
%  kalibrasi 1970-1999 (diambil dari hasil_koreksi_bias_hist_v2.mat),
%  lalu diterapkan pada data model CMIP6 2000-2014, dan hasilnya
%  dibandingkan dengan observasi ERA5-Land 2000-2014. Ini menguji
%  asumsi stasioneritas, apakah koreksi yang dilatih di satu periode
%  masih berlaku di periode lain (Teutschbein dan Seibert 2012,
%  differential split-sample test; Maraun 2016).
%
%  Efisiensi. Sisi kalibrasi tidak dihitung ulang. T_obs (ERA5
%  1970-1999) dan T_mod (CMIP6 1970-1999 hasil regrid) dibaca dari
%  file .mat kalibrasi. Skrip ini hanya membaca dan memproses data
%  2000-2014, jadi jauh lebih ringan dari pipeline penuh.
%
%  Keluaran disimpan dengan nama variabel yang sama seperti file
%  kalibrasi, sehingga skrip evaluasi (evaluasi_koreksi_bias.m) bisa
%  dijalankan tanpa diubah, cukup arahkan cfg.f_mat ke file ini.
%
%  Metode koreksi identik dengan pipeline, memakai QM.m, QDM.m,
%  fitter.m yang sama plus fallback empiris yang sama.
% ============================================================

clear; clc;

%% ---------------- KONFIGURASI ----------------
BASE = 'C:\Users\asraf\OneDrive\Dokumen\Akademik\Jurnal Fio\SKRIPSI FIORENZA (7)\SKRIPSI FIORENZA\Data Pulau Jawa';
cfg.f_cal_mat      = 'hasil_koreksi_bias_hist_v2.mat';
cfg.dir_cmip6_val  = fullfile(BASE, 'CMIP6', '2000-2014');
cfg.dir_era5_val   = fullfile(BASE, 'ERA5 NEW 2000-2014');
cfg.f_era5_val_mat = '';    % opsional, .mat berisi T_obs_val bulanan; kosongkan bila baca harian
cfg.f_out          = 'hasil_koreksi_bias_val_v2.mat';

cfg.thn_awal  = 2000;
cfg.thn_akhir = 2014;
cfg.nBulan    = (cfg.thn_akhir - cfg.thn_awal + 1) * 12;   % 180
cfg.win_qdm   = 120;
cfg.g         = 9.80665;

%% ---------------- MUAT SISI KALIBRASI ----------------
disp('Memuat sisi kalibrasi dari file .mat ...')
C = load(cfg.f_cal_mat, 'T_obs', 'T_mod', 'land_mask', 'elev_jawa', ...
    'idx_low', 'idx_high', 'prov_id', 'prov_names', 'lon_jawa', 'lat_jawa');
obs_cal   = C.T_obs;          % ERA5 1970-1999
mod_cal   = C.T_mod;          % CMIP6 1970-1999 regrid
land_mask = C.land_mask;
lon_jawa  = C.lon_jawa(:);
lat_jawa  = C.lat_jawa(:);
nLon = numel(lon_jawa);  nLat = numel(lat_jawa);
fprintf('Grid kalibrasi %d x %d, %d sel darat\n', nLon, nLat, sum(land_mask(:)))

%% ------------------------------------------------------------
%  1. OBSERVASI ERA5-LAND 2000-2014 (BULANAN)
% ------------------------------------------------------------
if ~isempty(cfg.f_era5_val_mat) && exist(cfg.f_era5_val_mat, 'file')
    disp('Tahap 1. Memuat ERA5 validasi bulanan dari .mat ...')
    S = load(cfg.f_era5_val_mat, 'T_obs_val');
    T_obs_val = S.T_obs_val;
    assert(isequal(size(T_obs_val), [nLon nLat cfg.nBulan]), ...
        'Ukuran T_obs_val tidak sesuai grid dan jumlah bulan')
else
    disp('Tahap 1. Ekstraksi ERA5-Land 2000-2014 harian ...')
    folders = dir(cfg.dir_era5_val);
    folders = folders([folders.isdir] & ~ismember({folders.name},{'.','..'}));

    T_daily = [];  t_daily = datetime.empty(0,1);  pertama = true;
    for y = 1:numel(folders)
        fls = dir(fullfile(cfg.dir_era5_val, folders(y).name, '*.nc'));
        for k = 1:numel(fls)
            f = fullfile(fls(k).folder, fls(k).name);
            [T, lonf, latf] = baca_nc_kanonik(f, 't2m');
            t = waktu_dari_nc(f, nama_dim_waktu(f, 't2m'));
            if pertama
                lon_idx = lonf >= min(lon_jawa)-1e-6 & lonf <= max(lon_jawa)+1e-6;
                lat_idx = latf >= min(lat_jawa)-1e-6 & latf <= max(lat_jawa)+1e-6;
                assert(sum(lon_idx)==nLon && sum(lat_idx)==nLat, ...
                    'Jendela ERA5 validasi tidak menghasilkan grid 91 x 30')
                assert(max(abs(lonf(lon_idx)-lon_jawa))<1e-4 && ...
                       max(abs(latf(lat_idx)-lat_jawa))<1e-4, ...
                    'Koordinat ERA5 validasi tidak sama dengan grid kalibrasi')
                pertama = false;
            end
            T_daily = cat(3, T_daily, T(lon_idx, lat_idx, :));
            t_daily = [t_daily; t(:)]; %#ok<AGROW>
        end
        fprintf('  selesai tahun %s\n', folders(y).name)
    end
    [t_daily, iord] = sort(t_daily);
    T_daily = T_daily(:,:,iord) - 273.15;

    disp('Agregasi harian ke bulanan ...')
    ym = (year(t_daily) - cfg.thn_awal) * 12 + month(t_daily);
    dalam = ym >= 1 & ym <= cfg.nBulan;
    T_obs_val = NaN(nLon, nLat, cfg.nBulan);
    for m = 1:cfg.nBulan
        sel = dalam & (ym == m);
        if any(sel), T_obs_val(:,:,m) = mean(T_daily(:,:,sel), 3, 'omitnan'); end
    end
    clear T_daily
end

%% ------------------------------------------------------------
%  2. MODEL CMIP6 2000-2014, REGRID BILINEAR KE GRID ERA5
% ------------------------------------------------------------
disp('Tahap 2. Ekstraksi dan regridding CMIP6 2000-2014 ...')
fls = dir(fullfile(cfg.dir_cmip6_val, '*.nc'));
assert(~isempty(fls), 'Tidak ada file .nc di %s', cfg.dir_cmip6_val)

T_raw = [];  t_mod = datetime.empty(0,1);  pertama = true;
for k = 1:numel(fls)
    f = fullfile(fls(k).folder, fls(k).name);
    [T, lon_c, lat_c] = baca_nc_kanonik(f, 'tas');
    t = waktu_dari_nc(f, nama_dim_waktu(f, 'tas'));
    if pertama
        mrg = 3 * abs(lon_c(2) - lon_c(1));
        src_i = lon_c >= min(lon_jawa)-mrg & lon_c <= max(lon_jawa)+mrg;
        src_j = lat_c >= min(lat_jawa)-mrg & lat_c <= max(lat_jawa)+mrg;
        lon_src = lon_c(src_i);  lat_src = lat_c(src_j);
        pertama = false;
    end
    T_raw = cat(3, T_raw, T(src_i, src_j, :));
    t_mod = [t_mod; t(:)]; %#ok<AGROW>
end
[t_mod, iord] = sort(t_mod);
T_raw = T_raw(:,:,iord) - 273.15;

ym_mod = (year(t_mod) - cfg.thn_awal) * 12 + month(t_mod);
sel = ym_mod >= 1 & ym_mod <= cfg.nBulan;
assert(numel(unique(ym_mod(sel))) == cfg.nBulan, ...
    'Bulan CMIP6 2000-2014 tidak lengkap 180')

T_mod_val = NaN(nLon, nLat, cfg.nBulan);
pos = find(sel);
for k = 1:numel(pos)
    V = squeeze(T_raw(:,:,pos(k)))';                 % (lat x lon)
    Tq = interp2(lon_src, lat_src, V, lon_jawa(:)', lat_jawa(:), 'linear');
    T_mod_val(:,:,ym_mod(pos(k))) = Tq';
end
clear T_raw
disp('Regridding selesai')

%% ------------------------------------------------------------
%  3. MASK VALIDITAS VALIDASI
% ------------------------------------------------------------
obs_lengkap = (sum(isfinite(T_obs_val), 3) == cfg.nBulan);
mod_lengkap = (sum(isfinite(T_mod_val), 3) == cfg.nBulan);
land_val = land_mask & obs_lengkap & mod_lengkap;

fprintf('Sel darat kalibrasi              = %d\n', sum(land_mask(:)))
fprintf('Sel valid pada 2000-2014         = %d\n', sum(land_val(:)))
fprintf('Sel gugur (ERA5 val tak lengkap) = %d\n', sum(land_mask(:) & ~obs_lengkap(:)))

%% ------------------------------------------------------------
%  4. TERAPKAN FUNGSI TRANSFER 1970-1999 KE 2000-2014
% ------------------------------------------------------------
disp('Tahap 4. Koreksi out-of-sample per grid ...')
T_QM_grid  = NaN(nLon, nLat, cfg.nBulan);
T_QDM_grid = NaN(nLon, nLat, cfg.nBulan);
stat = struct('n_sel',0,'qm_par',0,'qm_emp',0, ...
    'qdm_penuh',0,'qdm_tambal',0,'qdm_emp',0);

[si, sj] = find(land_val);
for c = 1:numel(si)
    i = si(c);  j = sj(c);
    o_cal = squeeze(obs_cal(i,j,:));     % 360, kalibrasi
    m_cal = squeeze(mod_cal(i,j,:));     % 360, kalibrasi
    x_val = squeeze(T_mod_val(i,j,:));   % 180, yang dikoreksi
    if any(~isfinite(m_cal)) || any(~isfinite(o_cal)) || any(~isfinite(x_val))
        continue
    end
    stat.n_sel = stat.n_sel + 1;

    % ---- QM parametrik dengan fallback empiris ----
    q = [];
    try, q = QM(o_cal, m_cal, x_val); catch, end
    if isempty(q) || any(~isfinite(q)) || std(q) < 1e-6
        q = eqm_empiris(o_cal, m_cal, x_val);
        stat.qm_emp = stat.qm_emp + 1;
    else
        stat.qm_par = stat.qm_par + 1;
    end
    T_QM_grid(i,j,:) = q;

    % ---- QDM parametrik, bulan gagal ditambal empiris ----
    try, cq = QDM(o_cal, m_cal, x_val); catch, cq = NaN(cfg.nBulan,1); end
    buruk = ~isfinite(cq);
    if any(buruk)
        cq(buruk) = eqdm_empiris(o_cal, m_cal, x_val, find(buruk), cfg.win_qdm);
        stat.qdm_tambal = stat.qdm_tambal + 1;
    else
        stat.qdm_penuh = stat.qdm_penuh + 1;
    end
    if std(cq,'omitnan') < 1e-6
        cq = eqdm_empiris(o_cal, m_cal, x_val, (1:cfg.nBulan)', cfg.win_qdm);
        stat.qdm_emp = stat.qdm_emp + 1;
    end
    T_QDM_grid(i,j,:) = cq;

    if mod(c,100)==0, fprintf('  %d dari %d sel\n', c, numel(si)); end
end

fprintf('\nRingkasan koreksi validasi\n')
fprintf('Sel diproses            = %d\n', stat.n_sel)
fprintf('QM parametrik / empiris = %d / %d\n', stat.qm_par, stat.qm_emp)
fprintf('QDM penuh / tambal / emp= %d / %d / %d\n', ...
    stat.qdm_penuh, stat.qdm_tambal, stat.qdm_emp)

% sanity cepat, bias regional out-of-sample
mask3 = repmat(land_val, 1, 1, cfg.nBulan);
b_raw = mean(T_mod_val(mask3) - T_obs_val(mask3), 'omitnan');
d_qm  = T_QM_grid  - T_obs_val;  b_qm  = mean(d_qm(mask3),  'omitnan');
d_qdm = T_QDM_grid - T_obs_val;  b_qdm = mean(d_qdm(mask3), 'omitnan');
fprintf('Bias regional 2000-2014  CMIP6 %.3f, QM %.3f, QDM %.3f (derajat C)\n', ...
    b_raw, b_qm, b_qdm)

%% ------------------------------------------------------------
%  5. SIMPAN DENGAN SKEMA YANG DIBACA SKRIP EVALUASI
% ------------------------------------------------------------
% nama variabel disamakan dengan file kalibrasi agar evaluasi_koreksi_bias.m
% bisa dijalankan tanpa diubah, cukup arahkan cfg.f_mat ke file ini
T_obs        = T_obs_val;             %#ok<NASGU>  observasi validasi
T_mod        = T_mod_val;             %#ok<NASGU>  model validasi (regrid)
land_mask    = land_val;              %#ok<NASGU>  mask valid validasi
elev_jawa    = C.elev_jawa;           %#ok<NASGU>
idx_low      = C.idx_low & land_val;  %#ok<NASGU>
idx_high     = C.idx_high & land_val; %#ok<NASGU>
prov_id      = C.prov_id;             %#ok<NASGU>
prov_names   = C.prov_names;          %#ok<NASGU>
time_bulanan = dateshift(datetime(cfg.thn_awal,1,1) + calmonths(0:cfg.nBulan-1), ...
    'end', 'month')';                 %#ok<NASGU>
diagnosa     = struct('stat', stat, 'land_val', land_val, ...
    'obs_lengkap', obs_lengkap, 'cfg', cfg); %#ok<NASGU>

save(cfg.f_out, 'T_QM_grid', 'T_QDM_grid', 'T_obs', 'T_mod', ...
    'land_mask', 'elev_jawa', 'idx_low', 'idx_high', ...
    'prov_id', 'prov_names', 'lon_jawa', 'lat_jawa', ...
    'time_bulanan', 'diagnosa', '-v7.3')
fprintf('\nHasil validasi tersimpan di %s\n', cfg.f_out)
fprintf('Jalankan evaluasi dengan mengarahkan cfg.f_mat ke %s\n', cfg.f_out)

%% ============================================================
%  FUNGSI LOKAL (identik dengan pipeline v2)
% ============================================================
function [data, lon, lat] = baca_nc_kanonik(fname, varname)
    info = ncinfo(fname, varname);
    dn = lower({info.Dimensions.Name});
    i_lon = find(ismember(dn, {'lon','longitude','x'}), 1);
    i_lat = find(ismember(dn, {'lat','latitude','y'}), 1);
    i_tim = find(ismember(dn, {'time','valid_time','t'}), 1);
    assert(~isempty(i_lon) && ~isempty(i_lat), ...
        'Dimensi lon lat tidak dikenali pada %s', fname)
    data = double(ncread(fname, varname));
    lon  = double(ncread(fname, info.Dimensions(i_lon).Name));
    lat  = double(ncread(fname, info.Dimensions(i_lat).Name));
    urut = [i_lon, i_lat, i_tim, setdiff(1:numel(dn), [i_lon i_lat i_tim])];
    data = permute(data, urut);
    if ndims(data) < 3
        data = reshape(data, size(data,1), size(data,2), []);
    end
    if lat(1) > lat(end), lat = flip(lat); data = flip(data, 2); end
    if lon(1) > lon(end), lon = flip(lon); data = flip(data, 1); end
end

function tname = nama_dim_waktu(fname, varname)
    info = ncinfo(fname, varname);
    dn = {info.Dimensions.Name};
    m = ismember(lower(dn), {'time','valid_time','t'});
    tname = dn{find(m, 1)};
end

function t = waktu_dari_nc(fname, tname)
    raw = double(ncread(fname, tname));
    units = ncreadatt(fname, tname, 'units');
    tok = regexp(units, '^\s*(\w+)\s+since\s+([\d\-]+)', 'tokens', 'once');
    ref = datetime(tok{2}, 'InputFormat', 'yyyy-M-d');
    switch lower(tok{1})
        case {'days','day'},        t = ref + days(raw);
        case {'hours','hour'},      t = ref + hours(raw);
        case {'minutes','minute'},  t = ref + minutes(raw);
        case {'seconds','second'},  t = ref + seconds(raw);
        otherwise, error('Satuan waktu %s tidak dikenali', tok{1})
    end
end

function xk = eqm_empiris(obs, mkal, x)
    p  = ((1:99)') / 100;
    qo = quantile(obs,  p);
    qm = quantile(mkal, p);
    [qm_u, iu] = unique(qm);
    qo_u = qo(iu);
    xk = interp1(qm_u, qo_u, x, 'linear', NaN);
    d_bwh = qo_u(1)   - qm_u(1);
    d_ats = qo_u(end) - qm_u(end);
    xk(x < qm_u(1))   = x(x < qm_u(1))   + d_bwh;
    xk(x > qm_u(end)) = x(x > qm_u(end)) + d_ats;
end

function c = eqdm_empiris(obs, mkal, x, idx_t, win)
    n = numel(x);
    c = NaN(numel(idx_t), 1);
    for k = 1:numel(idx_t)
        t = idx_t(k);
        a = max(1, t - win);
        b = min(n, t + win);
        xw = x(a:b);  xw = xw(isfinite(xw));
        tau = (sum(xw <= x(t)) - 0.5) / numel(xw);
        tau = min(max(tau, 1e-3), 1 - 1e-3);
        c(k) = quantile(obs, tau) + (x(t) - quantile(mkal, tau));
    end
end