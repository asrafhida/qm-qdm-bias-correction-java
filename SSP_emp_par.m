%% ============================================================
%  KOREKSI BIAS PROYEKSI CMIP6 PULAU JAWA
%  Skenario SSP2-4.5 dan SSP5-8.5, periode 2015 sampai 2100
%
%  PRINSIP
%  Fungsi transfer QM dan QDM dibangun HANYA dari periode kalibrasi
%  1970 sampai 1999, memakai ERA5-Land sebagai acuan dan CMIP6
%  historis sebagai model. Fungsi itu kemudian diterapkan pada luaran
%  CMIP6 periode proyeksi. Tidak ada observasi yang dipakai pada
%  periode proyeksi, sebab memang tidak tersedia.
%
%  DUA HAL PENTING YANG PERLU DIPAHAMI SEBELUM MENJALANKAN
%
%  1. Waktu jalan. QDM.m memanggil fitter di dalam loop untuk setiap
%     langkah waktu. Periode proyeksi berisi 1032 bulan, hampir tiga
%     kali panjang periode kalibrasi yang berisi 360 bulan. Karena
%     kalibrasi memakan sekitar 4 jam, satu skenario diperkirakan 9
%     sampai 12 jam dan dua skenario dapat melampaui 20 jam. Setel
%     cfg.metode_qdm ke 'empiris' untuk memakai QDM empiris berjendela
%     sama yang tidak melakukan pencocokan sebaran per langkah waktu,
%     sehingga jauh lebih cepat. Kuantil empiris justru merupakan cara
%     yang dipakai Cannon et al. (2015) dalam perumusan aslinya.
%
%  2. Pemotongan pada QM. Di dalam QM.m, peluang epst dijepit pada
%     1 dikurangi 1e-6 sebelum diinverskan. Untuk suhu proyeksi yang
%     melampaui rentang kalibrasi, terutama pada SSP5-8.5 akhir abad,
%     seluruh nilai ekstrem akan dipetakan ke satu nilai langit-langit
%     yang sama. Akibatnya sinyal pemanasan tertekan. Ini merupakan
%     wujud distorsi tren pada QM yang didokumentasikan Cannon et al.
%     (2015) dan Maraun (2016), dan justru menjadi alasan QDM
%     dirancang. Skrip ini menghitung diagnostik penjenuhan serta
%     perbandingan tren agar besarnya efek itu terukur, bukan sekadar
%     diasumsikan. Hasil QM pada periode proyeksi sebaiknya dilaporkan
%     sebagai pembanding metodologis, bukan sebagai proyeksi utama.
%
%  KELUARAN
%  Satu berkas .mat per skenario, berisi luaran mentah, QM, dan QDM
%  pada grid ERA5-Land, memakai skema nama variabel yang sama dengan
%  berkas kalibrasi sehingga skrip evaluasi dan peta dapat dipakai
%  ulang. Ditambah tabel tren dan tabel perubahan antarperiode.
%
%  Referensi
%  Cannon AJ, Sobie SR, Murdock TQ. 2015. J Climate 28:6938-6959.
%  Maraun D. 2016. Curr Clim Change Rep 2:211-220.
%  Muller WA et al. 2018. J Adv Model Earth Syst 10:1383-1413.
% ============================================================

clear; clc;

%% ---------------- KONFIGURASI ----------------
BASE = 'C:\Users\asraf\OneDrive\Dokumen\Akademik\Jurnal Fio\SKRIPSI FIORENZA (7)\SKRIPSI FIORENZA\Data Pulau Jawa';

cfg.f_cal_mat = 'hasil_koreksi_bias_hist_v2.mat';
cfg.skenario  = { ...
    'SSP2-4.5', fullfile(BASE,'CMIP6','ssp245'), 'hasil_proyeksi_ssp245_v2.mat'; ...
    'SSP5-8.5', fullfile(BASE,'CMIP6','ssp585'), 'hasil_proyeksi_ssp585_v2.mat'};

cfg.thn_awal    = 2015;
cfg.thn_akhir   = 2100;
cfg.metode_qdm  = 'parametrik';   % 'parametrik' memakai QDM.m, 'empiris' jauh lebih cepat
cfg.win_qdm     = 120;            % jendela bergerak, bulan, sama dengan QDM.m
cfg.f_xlsx      = 'ringkasan_proyeksi.xlsx';

% periode iris untuk pelaporan perubahan
cfg.periode = { ...
    'Dekat 2021-2040',  2021, 2040; ...
    'Menengah 2041-2060', 2041, 2060; ...
    'Jauh 2081-2100',   2081, 2100};

cfg.nBulan = (cfg.thn_akhir - cfg.thn_awal + 1) * 12;   % 1032

%% ---------------- MUAT SISI KALIBRASI ----------------
fprintf('Memuat sisi kalibrasi ...\n')
C = load(cfg.f_cal_mat, 'T_obs', 'T_mod', 'land_mask', 'elev_jawa', ...
    'idx_low', 'idx_high', 'prov_id', 'prov_names', 'lon_jawa', 'lat_jawa');
obs_cal = C.T_obs;
mod_cal = C.T_mod;
land    = logical(C.land_mask);
lon     = C.lon_jawa(:);
lat     = C.lat_jawa(:);
nLon = numel(lon);  nLat = numel(lat);
fprintf('Grid %d x %d, %d sel darat\n', nLon, nLat, sum(land(:)))

% acuan iklim masa kini dari observasi, dipakai sebagai dasar perubahan
base_obs = mean(obs_cal(repmat(land,1,1,size(obs_cal,3))), 'omitnan');
fprintf('Rata-rata ERA5-Land 1970-1999 di daratan Jawa = %.3f derajat C\n', base_obs)

if strcmpi(cfg.metode_qdm, 'parametrik')
    fprintf(['\nPERINGATAN. cfg.metode_qdm = parametrik. QDM.m mencocokkan\n' ...
        'sebaran untuk tiap langkah waktu, sehingga %d bulan kali %d sel\n' ...
        'berarti sekitar %.1f juta pemanggilan fitter per skenario.\n' ...
        'Perkiraan waktu 9 sampai 12 jam per skenario. Setel ke empiris\n' ...
        'bila ingin jauh lebih cepat.\n\n'], cfg.nBulan, sum(land(:)), ...
        cfg.nBulan*sum(land(:))/1e6)
end

%% ---------------- PROSES TIAP SKENARIO ----------------
nSk = size(cfg.skenario, 1);
Tren = table();  Ubah = table();

for s = 1:nSk
    nama_sk = cfg.skenario{s,1};
    dir_sk  = cfg.skenario{s,2};
    f_out   = cfg.skenario{s,3};
    fprintf('\n================ %s ================\n', nama_sk)

    %% 1. Ekstraksi dan regridding
    fprintf('Ekstraksi dan regridding CMIP6 ...\n')
    fls = dir(fullfile(dir_sk, '*.nc'));
    assert(~isempty(fls), 'Tidak ada berkas .nc di %s', dir_sk)

    T_raw = [];  t_mod = datetime.empty(0,1);  pertama = true;
    for k = 1:numel(fls)
        f = fullfile(fls(k).folder, fls(k).name);
        [T, lon_c, lat_c] = baca_nc_kanonik(f, 'tas');
        t = waktu_dari_nc(f, nama_dim_waktu(f, 'tas'));
        if pertama
            mrg = 3 * abs(lon_c(2) - lon_c(1));
            si = lon_c >= min(lon)-mrg & lon_c <= max(lon)+mrg;
            sj = lat_c >= min(lat)-mrg & lat_c <= max(lat)+mrg;
            lon_src = lon_c(si);  lat_src = lat_c(sj);
            pertama = false;
        end
        T_raw = cat(3, T_raw, T(si, sj, :));
        t_mod = [t_mod; t(:)]; %#ok<AGROW>
    end
    [t_mod, iord] = sort(t_mod);
    T_raw = T_raw(:,:,iord) - 273.15;

    fprintf('Rentang berkas %d sampai %d\n', year(t_mod(1)), year(t_mod(end)))
    ym = (year(t_mod) - cfg.thn_awal) * 12 + month(t_mod);
    sel = ym >= 1 & ym <= cfg.nBulan;
    n_unik = numel(unique(ym(sel)));
    if n_unik < cfg.nBulan
        warning('Hanya %d dari %d bulan tersedia, sisanya akan NaN', ...
            n_unik, cfg.nBulan)
    end

    T_proy = NaN(nLon, nLat, cfg.nBulan);
    pos = find(sel);
    for k = 1:numel(pos)
        V  = squeeze(T_raw(:,:,pos(k)))';
        Tq = interp2(lon_src, lat_src, V, lon(:)', lat(:), 'linear');
        T_proy(:,:,ym(pos(k))) = Tq';
    end
    clear T_raw
    fprintf('Regridding selesai\n')

    %% 2. Mask sel yang dapat dikoreksi
    proy_lengkap = (sum(isfinite(T_proy), 3) == cfg.nBulan);
    land_s = land & proy_lengkap;
    fprintf('Sel dapat dikoreksi = %d dari %d sel darat\n', ...
        sum(land_s(:)), sum(land(:)))

    %% 3. Terapkan fungsi transfer kalibrasi
    fprintf('Koreksi bias proyeksi, metode QDM %s ...\n', cfg.metode_qdm)
    T_QM_grid  = NaN(nLon, nLat, cfg.nBulan);
    T_QDM_grid = NaN(nLon, nLat, cfg.nBulan);
    stat = struct('n',0,'qm_par',0,'qm_emp',0,'qdm_par',0,'qdm_emp',0, ...
        'qm_jenuh',0);

    [si2, sj2] = find(land_s);
    t0 = tic;
    for c = 1:numel(si2)
        i = si2(c);  j = sj2(c);
        o_cal = squeeze(obs_cal(i,j,:));
        m_cal = squeeze(mod_cal(i,j,:));
        x_pro = squeeze(T_proy(i,j,:));
        if any(~isfinite(o_cal)) || any(~isfinite(m_cal)), continue, end
        stat.n = stat.n + 1;

        % ---- QM ----
        q = [];
        try, q = QM(o_cal, m_cal, x_pro); catch, end
        if isempty(q) || any(~isfinite(q)) || std(q) < 1e-6
            q = eqm_empiris(o_cal, m_cal, x_pro);
            stat.qm_emp = stat.qm_emp + 1;
        else
            stat.qm_par = stat.qm_par + 1;
        end
        T_QM_grid(i,j,:) = q;

        % diagnostik penjenuhan, berapa bulan menempel pada nilai tertinggi
        tol = 1e-4;
        stat.qm_jenuh = stat.qm_jenuh + sum(abs(q - max(q)) < tol) - 1;

        % ---- QDM ----
        if strcmpi(cfg.metode_qdm, 'parametrik')
            try, cq = QDM(o_cal, m_cal, x_pro); catch, cq = NaN(cfg.nBulan,1); end
            buruk = ~isfinite(cq);
            if any(buruk)
                cq(buruk) = eqdm_empiris(o_cal, m_cal, x_pro, find(buruk), cfg.win_qdm);
            end
            stat.qdm_par = stat.qdm_par + 1;
        else
            cq = eqdm_empiris(o_cal, m_cal, x_pro, (1:cfg.nBulan)', cfg.win_qdm);
            stat.qdm_emp = stat.qdm_emp + 1;
        end
        T_QDM_grid(i,j,:) = cq;

        if mod(c, 10) == 0
            el = toc(t0);
            fprintf('  %d dari %d sel, %.1f menit berjalan, perkiraan sisa %.1f menit\n', ...
                c, numel(si2), el/60, el/60*(numel(si2)-c)/c);
        end
    end
    fprintf('Selesai dalam %.1f menit\n', toc(t0)/60)
    fprintf('QM parametrik %d, QM empiris %d, QDM %d\n', ...
        stat.qm_par, stat.qm_emp, stat.qdm_par + stat.qdm_emp)
    fprintf('Rata-rata bulan menempel nilai puncak QM per sel = %.1f dari %d\n', ...
        stat.qm_jenuh / max(1,stat.n), cfg.nBulan)

    %% 4. Diagnostik tren dan perubahan antarperiode
    thn = repelem((cfg.thn_awal:cfg.thn_akhir)', 12);
    r_raw = regmean(T_proy,     land_s, cfg.nBulan);
    r_qm  = regmean(T_QM_grid,  land_s, cfg.nBulan);
    r_qdm = regmean(T_QDM_grid, land_s, cfg.nBulan);

    [tr_raw, a_raw] = tren_dekade(r_raw, thn);
    [tr_qm,  a_qm ] = tren_dekade(r_qm,  thn);
    [tr_qdm, a_qdm] = tren_dekade(r_qdm, thn);

    fprintf('\n--- Tren pemanasan regional, derajat C per dekade ---\n')
    fprintf('CMIP6 mentah %.4f, QM %.4f, QDM %.4f\n', tr_raw, tr_qm, tr_qdm)
    fprintf('QM mempertahankan %.1f persen tren mentah, QDM %.1f persen\n', ...
        100*tr_qm/tr_raw, 100*tr_qdm/tr_raw)

    Tren = [Tren; table(string(nama_sk), tr_raw, tr_qm, tr_qdm, ...
        100*tr_qm/tr_raw, 100*tr_qdm/tr_raw, ...
        'VariableNames', {'Skenario','Tren_CMIP6','Tren_QM','Tren_QDM', ...
        'Persen_QM','Persen_QDM'})]; %#ok<AGROW>

    fprintf('\n--- Suhu rata-rata daratan Jawa per periode ---\n')
    fprintf('%-20s %10s %10s %10s %12s\n', 'Periode', 'CMIP6', 'QM', 'QDM', 'dQDM_thd_ERA5')
    for pp = 1:size(cfg.periode,1)
        m = thn >= cfg.periode{pp,2} & thn <= cfg.periode{pp,3};
        v_raw = mean(r_raw(m), 'omitnan');
        v_qm  = mean(r_qm(m),  'omitnan');
        v_qdm = mean(r_qdm(m), 'omitnan');
        fprintf('%-20s %10.3f %10.3f %10.3f %12.3f\n', cfg.periode{pp,1}, ...
            v_raw, v_qm, v_qdm, v_qdm - base_obs)
        Ubah = [Ubah; table(string(nama_sk), string(cfg.periode{pp,1}), ...
            v_raw, v_qm, v_qdm, v_qdm - base_obs, v_qm - base_obs, ...
            'VariableNames', {'Skenario','Periode','CMIP6','QM','QDM', ...
            'dQDM','dQM'})]; %#ok<AGROW>
    end

    % tren per sel, untuk dipetakan kemudian
    tren_sel = NaN(nLon, nLat, 3);
    for c = 1:numel(si2)
        i = si2(c);  j = sj2(c);
        tren_sel(i,j,1) = tren_dekade(squeeze(T_proy(i,j,:)),     thn);
        tren_sel(i,j,2) = tren_dekade(squeeze(T_QM_grid(i,j,:)),  thn);
        tren_sel(i,j,3) = tren_dekade(squeeze(T_QDM_grid(i,j,:)), thn);
    end

    %% 5. Simpan dengan skema yang sama seperti berkas kalibrasi
    T_mod        = T_proy;               %#ok<NASGU>
    land_mask    = land_s;               %#ok<NASGU>
    elev_jawa    = C.elev_jawa;          %#ok<NASGU>
    idx_low      = C.idx_low  & land_s;  %#ok<NASGU>
    idx_high     = C.idx_high & land_s;  %#ok<NASGU>
    prov_id      = C.prov_id;            %#ok<NASGU>
    prov_names   = C.prov_names;         %#ok<NASGU>
    lon_jawa     = lon;                  %#ok<NASGU>
    lat_jawa     = lat;                  %#ok<NASGU>
    time_bulanan = dateshift(datetime(cfg.thn_awal,1,1) + ...
        calmonths(0:cfg.nBulan-1), 'end', 'month')';   %#ok<NASGU>
    diagnosa     = struct('stat', stat, 'skenario', nama_sk, ...
        'tren_regional', [tr_raw tr_qm tr_qdm], 'tren_sel', tren_sel, ...
        'base_obs', base_obs, 'cfg', cfg);             %#ok<NASGU>

    save(f_out, 'T_QM_grid', 'T_QDM_grid', 'T_mod', ...
        'land_mask', 'elev_jawa', 'idx_low', 'idx_high', ...
        'prov_id', 'prov_names', 'lon_jawa', 'lat_jawa', ...
        'time_bulanan', 'diagnosa', '-v7.3')
    fprintf('\nTersimpan di %s\n', f_out)

    clear T_proy T_QM_grid T_QDM_grid T_mod
end

%% ---------------- SIMPAN RINGKASAN ----------------
writetable(Tren, cfg.f_xlsx, 'Sheet','Tren')
writetable(Ubah, cfg.f_xlsx, 'Sheet','Perubahan')
fprintf('\nRingkasan tersimpan di %s\n', cfg.f_xlsx)
fprintf(['\nCatatan pembacaan. Bila persentase tren QM jauh di bawah 100\n' ...
    'sedangkan QDM mendekati 100, itu adalah bukti kuantitatif distorsi\n' ...
    'tren pada QM dan pemeliharaan tren pada QDM, sesuai Cannon et al.\n' ...
    '(2015). Efeknya diperkirakan lebih besar pada SSP5-8.5 karena\n' ...
    'suhu proyeksinya lebih jauh melampaui rentang kalibrasi.\n'])

%% ============================================================
%  FUNGSI LOKAL
% ============================================================
function [tr, a] = tren_dekade(v, thn)
% Tren linear atas rata-rata tahunan, dinyatakan dalam derajat C per
% dekade. Rata-rata tahunan dipakai lebih dahulu agar siklus musiman
% tidak membebani pendugaan kemiringan.
    u = unique(thn);
    ya = arrayfun(@(y) mean(v(thn == y), 'omitnan'), u);
    ok = isfinite(ya);
    if nnz(ok) < 10, tr = NaN;  a = NaN;  return, end
    p = polyfit(u(ok), ya(ok), 1);
    tr = p(1) * 10;
    a  = p(2);
end

function r = regmean(X, land, nT)
    m = repmat(land, 1, 1, nT);
    X(~m) = NaN;
    r = squeeze(mean(mean(X, 1, 'omitnan'), 2, 'omitnan'));
end

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
    if lat(1) > lat(end), lat = flip(lat);  data = flip(data, 2); end
    if lon(1) > lon(end), lon = flip(lon);  data = flip(data, 1); end
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
        case {'days','day'},       t = ref + days(raw);
        case {'hours','hour'},     t = ref + hours(raw);
        case {'minutes','minute'}, t = ref + minutes(raw);
        case {'seconds','second'}, t = ref + seconds(raw);
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
% QDM empiris aditif mengikuti struktur persamaan Cannon et al. (2015).
% Posisi kuantil dihitung pada jendela bergerak deret proyeksi, lalu
% nilai observasi pada kuantil itu dijumlahkan dengan selisih terhadap
% model kalibrasi. Suku selisih inilah yang memelihara sinyal
% perubahan sehingga tren proyeksi tidak tertekan.
    n = numel(x);
    c = NaN(numel(idx_t), 1);
    for k = 1:numel(idx_t)
        t = idx_t(k);
        a = max(1, t - win);
        b = min(n, t + win);
        xw = x(a:b);  xw = xw(isfinite(xw));
        if isempty(xw), continue, end
        tau = (sum(xw <= x(t)) - 0.5) / numel(xw);
        tau = min(max(tau, 1e-3), 1 - 1e-3);
        c(k) = quantile(obs, tau) + (x(t) - quantile(mkal, tau));
    end
end