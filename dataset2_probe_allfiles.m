%% Dataset 2 Analysis 
% Capstone II - Team E(Chidna)
% JiaJia Fu 
% 4/20/26

clear;
clc;
close all;

%% --- PARAMETERS ---
f0              = 0.5e6;
harmonics       = 6:9;           % matches single-frame: 6:9
n_harmonics     = length(harmonics);
uh_orders       = [6.5 7.5 8.5];
bin_range       = 1000;
threshold_N     = 1;

%% --- PATH SETUP ---
function_path = '/Users/jiasquared/Desktop/CODING/capstone_team_echidna/04_Projects/';
addpath(function_path);

probe_path   = '/Users/jiasquared/Desktop/CODING/capstone_team_echidna/1_and_2_2026-01-16_single_tube_sweep/probe';

probe_files  = dir(fullfile(probe_path, '*.pacq*'));
target_files = probe_files;
ambient_file_idx = find(contains({probe_files.name}, 'ambient'), 1);

if isempty(ambient_file_idx)
    error('No ambient .pacq file found');
end

disp(['Ambient file detected: ', probe_files(ambient_file_idx).name]);

%% --- PREALLOCATE RESULTS ---
N = length(target_files);

Results(N) = struct( ...
    'filename',                [], ...
    'frame',                   [], ...
    'harmonic_amplitudes',     [], ...   % raw PSD values at H freqs (matches WelchSpec(harmonic_indices))
    'target_harm_dB',          [], ...   % 10*log10 of above
    'ultraharmonic_amplitudes',[], ...   % raw PSD values at UH freqs
    'target_uh_dB',            [], ...   % 10*log10 of above
    'harmonic_indices',        [], ...
    'ultraharmonic_indices',   [], ...
    'exceeds_harm',            [], ...
    'exceeds_uh',              [], ...
    'harmonic_threshold_dB',   [], ...
    'uh_threshold_dB',         [], ...
    'target_bb_std',           [], ...
    'ambient_bb_std',          [], ...
    'safety_1s',               [], ...
    'safety_2s',               [], ...
    'safety_3s',               []);

%% --- MAIN LOOP: all frames ---
N_files = numel(target_files);
Results = [];
 
% --- Load ambient data once ---
ambient_data = extract_probe_data(probe_path, ambient_file_idx);
a_fs     = ambient_data.fs;
N_ambient = size(ambient_data.rf_data, 3);
 
fprintf('Loaded ambient file: %s (%d frames)\n', ...
    target_files(ambient_file_idx).name, N_ambient);
 
% --- Ambient diagnostics ---
for f = 1:N_ambient
    frame_amb     = double(ambient_data.rf_data(:, 64, f));
    mean_pressure = mean(frame_amb);
    std_pressure  = std(frame_amb);
    max_pressure  = max(frame_amb);
    min_pressure  = min(frame_amb);
 
    [pxx_a, f_vec] = pwelch(frame_amb, [], [], [], a_fs);
    df_amb = f_vec(2) - f_vec(1);
 
    bb_lo_hz = 6.5 * f0 + bin_range * df_amb;
    bb_hi_hz = 7.0 * f0 - bin_range * df_amb;
    bb_lo_n  = max(min(bb_lo_hz, bb_hi_hz), 1e-4) / (a_fs/2);
    bb_hi_n  = min(max(bb_lo_hz, bb_hi_hz), a_fs/2 - 1) / (a_fs/2);
    bb_hi_n  = min(bb_hi_n, 1-1e-4);
    [b_bp, a_bp] = butter(4, [bb_lo_n, bb_hi_n], 'bandpass');
    bb_std = std(filtfilt(b_bp, a_bp, frame_amb));
 
    fprintf('Ambient Frame %d: Mean=%.3f Std=%.3f Max=%.3f Min=%.3f BB_std=%.3f\n', ...
        f, mean_pressure, std_pressure, max_pressure, min_pressure, bb_std);
end
 
%% --- Per-file, per-frame loop ---
for i = 1:N_files
    if i == ambient_file_idx
        continue
    end
 
    fprintf('\n[%d/%d] Processing: %s\n', i, N_files, target_files(i).name);
 
    clean_data = extract_probe_data(probe_path, i);
    fs_i   = clean_data.fs;
    N_frames = size(clean_data.rf_data, 3);
 
    % --- Parse concentration label from filename ---
    % Expects filenames like PCI_1000x_01.pacq → concentration = '1000x_01'
    raw_name = target_files(i).name;
    conc_token = regexp(raw_name, 'PCI_(.+?)\.', 'tokens', 'once');
    if ~isempty(conc_token)
        concentration = conc_token{1};
    else
        concentration = raw_name;  % fallback
    end
 
    harmonic_freqs = f0 * harmonics;  % length = n_harmonics (3)
 
    for f = 1:N_frames
 
        frame_i         = double(clean_data.rf_data(:, 64, f));
        ambient_frame_i = double(ambient_data.rf_data(:, 64, f));
 
        % --- Pressure from frame index ---
        % Frame 1 = 0.025 MPa, increments of 0.005 MPa
        pressure_MPa = 0.020 + f * 0.005;
 
        % --- PSDs ---
        [pxx_i,   f_i]   = pwelch(frame_i,         [], [], [], fs_i);
        [pxx_amb, f_amb]  = pwelch(ambient_frame_i, [], [], [], a_fs);
        df_i   = f_i(2)   - f_i(1);
        df_amb = f_amb(2)  - f_amb(1);
        search_range_i   = bin_range * df_i;
        search_range_amb = bin_range * df_amb;
 
        % ---- Ambient thresholds ----------------------------------------
        % Ultraharmonics
        amb_uh_peak_dB = NaN(1, length(uh_orders));
        amb_uh_std_dB  = NaN(1, length(uh_orders));
        for k = 1:length(uh_orders)
            f_target   = uh_orders(k) * f0;
            idx_w = find(f_amb >= f_target - search_range_amb & ...
                         f_amb <= f_target + search_range_amb);
            if ~isempty(idx_w)
                amb_uh_peak_dB(k) = 10*log10(max(pxx_amb(idx_w)));
                amb_uh_std_dB(k)  = std(10*log10(pxx_amb(idx_w) + eps));
            end
        end
        uh_threshold_dB = amb_uh_peak_dB + threshold_N * amb_uh_std_dB;
 
        % Harmonics noise floor
        amb_noise_mask = true(size(f_amb));
        for k = 1:n_harmonics
            [~, h_idx] = min(abs(f_amb - harmonic_freqs(k)));
            lo = max(1, h_idx - bin_range);
            hi = min(length(f_amb), h_idx + bin_range);
            amb_noise_mask(lo:hi) = false;
        end
        for k = 1:length(uh_orders)
            [~, u_idx] = min(abs(f_amb - uh_orders(k)*f0));
            lo = max(1, u_idx - bin_range);
            hi = min(length(f_amb), u_idx + bin_range);
            amb_noise_mask(lo:hi) = false;
        end
        amb_noise_vals_dB     = 10*log10(pxx_amb(amb_noise_mask));
        amb_noise_mean_dB     = mean(amb_noise_vals_dB);
        amb_noise_std_dB      = std(amb_noise_vals_dB);
        harmonic_threshold_dB = amb_noise_mean_dB + threshold_N * amb_noise_std_dB;
        % scalar threshold — same for all harmonics (matches PCD broadband approach)
 
        % ---- Broadband filter ------------------------------------------
        bb_lo_hz = 6.5 * f0 + bin_range * df_amb;
        bb_hi_hz = 7.0 * f0 - bin_range * df_amb;
        bb_lo_i  = max(min(bb_lo_hz, bb_hi_hz), 1e-4) / (fs_i/2);
        bb_hi_i  = min(max(bb_lo_hz, bb_hi_hz), fs_i/2 - 1) / (fs_i/2);
        bb_hi_i  = min(bb_hi_i, 1-1e-4);
        [b_i, a_i] = butter(4, [bb_lo_i, bb_hi_i], 'bandpass');
 
        ambient_bb_std  = std(filtfilt(b_i, a_i, ambient_frame_i));
        target_bb_std_i = std(filtfilt(b_i, a_i, frame_i));
 
        % ---- Harmonics (3): H1, H2, H3 --------------------------------
        % Find index of each harmonic in f_i
        h_idx_vec = arrayfun(@(hf) find(abs(f_i - hf) == min(abs(f_i - hf)), 1), ...
                             harmonic_freqs);
        h_amp_lin = pxx_i(h_idx_vec);   % linear power at peak bin
        h_dB      = 10*log10(h_amp_lin);
 
        % ---- Ultraharmonics (3): U1, U2, U3 ---------------------------
        uh_idx_vec = zeros(1, numel(uh_orders));
        uh_amp_lin = zeros(1, numel(uh_orders));
        for k = 1:numel(uh_orders)
            f_target = uh_orders(k) * f0;
            idx_w = find(f_i >= f_target - search_range_i & ...
                         f_i <= f_target + search_range_i);
            if ~isempty(idx_w)
                [max_val, max_loc] = max(pxx_i(idx_w));
                uh_idx_vec(k) = idx_w(max_loc);
                uh_amp_lin(k) = max_val;
            end
        end
        uh_dB = 10*log10(uh_amp_lin);
 
        % ---- Threshold flags -------------------------------------------
        exceeds_harm = h_dB  > harmonic_threshold_dB;   % 1×3 logical
        exceeds_uh   = uh_dB > uh_threshold_dB;         % 1×3 logical
 
        % ---- Safety flags ----------------------------------------------
        safety_1s = target_bb_std_i > 2 * ambient_bb_std;
        safety_2s = target_bb_std_i > 3 * ambient_bb_std;
        safety_3s = target_bb_std_i > 4 * ambient_bb_std;
 
        % ---- Console output --------------------------------------------
        fprintf('  Frame %d/%d | Pressure %.4f MPa\n', f, N_frames, pressure_MPa);
        fprintf('    H_dB:  [%.2f  %.2f  %.2f]\n',  h_dB(1),  h_dB(2),  h_dB(3));
        fprintf('    UH_dB: [%.2f  %.2f  %.2f]\n', uh_dB(1), uh_dB(2), uh_dB(3));
        fprintf('    Exceeds harm: [%d %d %d]  uh: [%d %d %d]\n', ...
            exceeds_harm(1), exceeds_harm(2), exceeds_harm(3), ...
            exceeds_uh(1),   exceeds_uh(2),   exceeds_uh(3));
        fprintf('    BB std: %.4f  Ambient: %.4f  Safety: 1s=%d 2s=%d 3s=%d\n', ...
            target_bb_std_i, ambient_bb_std, safety_1s, safety_2s, safety_3s);
 
        % ---- Store as flat struct (one scalar/string per field) --------
        % Pre-define ALL fields with defaults so struct shape is identical
        % every iteration — prevents vertcat field-mismatch errors.
        R = struct( ...
            'concentration',         '', ...
            'filename',              '', ...
            'frame',                 NaN, ...
            'pressure_MPa',          NaN, ...
            'H1_idx',                NaN, ...
            'H2_idx',                NaN, ...
            'H3_idx',                NaN, ...
            'H1_dB',                 NaN, ...
            'H2_dB',                 NaN, ...
            'H3_dB',                 NaN, ...
            'U1_idx',                NaN, ...
            'U2_idx',                NaN, ...
            'U3_idx',                NaN, ...
            'U1_dB',                 NaN, ...
            'U2_dB',                 NaN, ...
            'U3_dB',                 NaN, ...
            'target_bb_std',         NaN, ...
            'ambient_bb_std',        NaN, ...
            'harmonic_threshold_dB', NaN, ...
            'U1_threshold_dB',       NaN, ...
            'U2_threshold_dB',       NaN, ...
            'U3_threshold_dB',       NaN, ...
            'exceeds_H1',            false, ...
            'exceeds_H2',            false, ...
            'exceeds_H3',            false, ...
            'exceeds_U1',            false, ...
            'exceeds_U2',            false, ...
            'exceeds_U3',            false, ...
            'safety_1s',             false, ...
            'safety_2s',             false, ...
            'safety_3s',             false  ...
        );
 
        R.concentration     = concentration;
        R.filename          = raw_name;
        R.frame             = f;
        R.pressure_MPa      = pressure_MPa;
 
        % Harmonic indices
        R.H1_idx            = h_idx_vec(1);
        R.H2_idx            = h_idx_vec(2);
        R.H3_idx            = h_idx_vec(3);
 
        % Harmonic dB
        R.H1_dB             = h_dB(1);
        R.H2_dB             = h_dB(2);
        R.H3_dB             = h_dB(3);
 
        % Ultraharmonic indices
        R.U1_idx            = uh_idx_vec(1);
        R.U2_idx            = uh_idx_vec(2);
        R.U3_idx            = uh_idx_vec(3);
 
        % Ultraharmonic dB
        R.U1_dB             = uh_dB(1);
        R.U2_dB             = uh_dB(2);
        R.U3_dB             = uh_dB(3);
 
        % Broadband
        R.target_bb_std     = target_bb_std_i;
        R.ambient_bb_std    = ambient_bb_std;
 
        % Threshold values used
        R.harmonic_threshold_dB = harmonic_threshold_dB;   % scalar
        R.U1_threshold_dB   = uh_threshold_dB(1);
        R.U2_threshold_dB   = uh_threshold_dB(2);
        R.U3_threshold_dB   = uh_threshold_dB(3);
 
        % Threshold flags
        R.exceeds_H1        = exceeds_harm(1);
        R.exceeds_H2        = exceeds_harm(2);
        R.exceeds_H3        = exceeds_harm(3);
        R.exceeds_U1        = exceeds_uh(1);
        R.exceeds_U2        = exceeds_uh(2);
        R.exceeds_U3        = exceeds_uh(3);
 
        % Safety flags
        R.safety_1s         = safety_1s;
        R.safety_2s         = safety_2s;
        R.safety_3s         = safety_3s;
 
        Results = [Results; R];
    end
end
 
%% --- Save probe CSV ---
T_probe = struct2table(Results);
writetable(T_probe, 'Probe_Features.csv');
fprintf('\nProbe CSV saved: %d rows.\n', height(T_probe));
 

%% Functions

function out = extract_probe_data(probe_path, i)
    data_files = dir(fullfile(probe_path, '*.pacq*'));
    file_path  = string(fullfile(data_files(i).folder, data_files(i).name));
    acq        = imaging.common.acquisition_reader(file_path);
    disp(acq.data_params)
    disp(acq)

    N_t      = acq.rfdata_dimensions(1);
    N_ch     = acq.rfdata_dimensions(2);
    N_frames = acq.number_of_frames;
    fs       = acq.data_params.fs;

    rf_data  = zeros(N_t, N_ch, N_frames, 'double');
    for k = 1:N_frames
        rf_data(:, :, k) = double(acq.data{k});
    end

    out.rf_data  = rf_data;
    out.t        = (0:N_t-1).' / fs;
    out.frame_t  = acq.elapsed_timevector(:);
    out.fs       = fs;
end

function features = get_features(f, pxx, f0)
    pxx_db = 10*log10(pxx);
    df     = f(2) - f(1);

    harmonics      = 6:12;
    ultraharmonics = 6.5:1:12.5;
    search_window  = round(0.08e6 / df);

    features        = struct();
    harm_freq_vals  = zeros(1, length(harmonics));
    ultra_freq_vals = zeros(1, length(ultraharmonics));
end fprintf('Loaded ambient file: %s (%d frames)\n', ...
    target_files(ambient_file_idx).name, N_ambient);
 
% --- Ambient diagnostics ---
for f = 1:N_ambient
    frame_amb     = double(ambient_data.rf_data(:, 64, f));
    mean_pressure = mean(frame_amb);
    std_pressure  = std(frame_amb);
    max_pressure  = max(frame_amb);
    min_pressure  = min(frame_amb);
 
    [pxx_a, f_vec] = pwelch(frame_amb, [], [], [], a_fs);
    df_amb = f_vec(2) - f_vec(1);
 
    bb_lo_hz = 6.5 * f0 + bin_range * df_amb;
    bb_hi_hz = 7.0 * f0 - bin_range * df_amb;
    bb_lo_n  = max(min(bb_lo_hz, bb_hi_hz), 1e-4) / (a_fs/2);
    bb_hi_n  = min(max(bb_lo_hz, bb_hi_hz), a_fs/2 - 1) / (a_fs/2);
    bb_hi_n  = min(bb_hi_n, 1-1e-4);
    [b_bp, a_bp] = butter(4, [bb_lo_n, bb_hi_n], 'bandpass');
    bb_std = std(filtfilt(b_bp, a_bp, frame_amb));
 
    fprintf('Ambient Frame %d: Mean=%.3f Std=%.3f Max=%.3f Min=%.3f BB_std=%.3f\n', ...
        f, mean_pressure, std_pressure, max_pressure, min_pressure, bb_std);
end
 
%% --- Per-file, per-frame loop ---
for i = 1:N_files
    if i == ambient_file_idx
        continue
    end
 
    fprintf('\n[%d/%d] Processing: %s\n', i, N_files, target_files(i).name);
 
    clean_data = extract_probe_data(probe_path, i);
    fs_i   = clean_data.fs;
    N_frames = size(clean_data.rf_data, 3);
 
    % --- Parse concentration label from filename ---
    % Expects filenames like PCI_1000x_01.pacq → concentration = '1000x_01'
    raw_name = target_files(i).name;
    conc_token = regexp(raw_name, 'PCI_(.+?)\.', 'tokens', 'once');
    if ~isempty(conc_token)
        concentration = conc_token{1};
    else
        concentration = raw_name;  % fallback
    end
 
    harmonic_freqs = f0 * harmonics;  % length = n_harmonics (3)
 
    for f = 1:N_frames
 
        frame_i         = double(clean_data.rf_data(:, 64, f));
        ambient_frame_i = double(ambient_data.rf_data(:, 64, f));
 
        % --- Pressure from frame index ---
        % Frame 1 = 0.025 MPa, increments of 0.005 MPa
        pressure_MPa = 0.020 + f * 0.005;
 
        % --- PSDs ---
        [pxx_i,   f_i]   = pwelch(frame_i,         [], [], [], fs_i);
        [pxx_amb, f_amb]  = pwelch(ambient_frame_i, [], [], [], a_fs);
        df_i   = f_i(2)   - f_i(1);
        df_amb = f_amb(2)  - f_amb(1);
        search_range_i   = bin_range * df_i;
        search_range_amb = bin_range * df_amb;
 
        % ---- Ambient thresholds ----------------------------------------
        % Ultraharmonics
        amb_uh_peak_dB = NaN(1, length(uh_orders));
        amb_uh_std_dB  = NaN(1, length(uh_orders));
        for k = 1:length(uh_orders)
            f_target   = uh_orders(k) * f0;
            idx_w = find(f_amb >= f_target - search_range_amb & ...
                         f_amb <= f_target + search_range_amb);
            if ~isempty(idx_w)
                amb_uh_peak_dB(k) = 10*log10(max(pxx_amb(idx_w)));
                amb_uh_std_dB(k)  = std(10*log10(pxx_amb(idx_w) + eps));
            end
        end
        uh_threshold_dB = amb_uh_peak_dB + threshold_N * amb_uh_std_dB;
 
        % Harmonics noise floor
        amb_noise_mask = true(size(f_amb));
        for k = 1:n_harmonics
            [~, h_idx] = min(abs(f_amb - harmonic_freqs(k)));
            lo = max(1, h_idx - bin_range);
            hi = min(length(f_amb), h_idx + bin_range);
            amb_noise_mask(lo:hi) = false;
        end
        for k = 1:length(uh_orders)
            [~, u_idx] = min(abs(f_amb - uh_orders(k)*f0));
            lo = max(1, u_idx - bin_range);
            hi = min(length(f_amb), u_idx + bin_range);
            amb_noise_mask(lo:hi) = false;
        end
        amb_noise_vals_dB     = 10*log10(pxx_amb(amb_noise_mask));
        amb_noise_mean_dB     = mean(amb_noise_vals_dB);
        amb_noise_std_dB      = std(amb_noise_vals_dB);
        harmonic_threshold_dB = amb_noise_mean_dB + threshold_N * amb_noise_std_dB;
        % scalar threshold — same for all harmonics (matches PCD broadband approach)
 
        % ---- Broadband filter ------------------------------------------
        bb_lo_hz = 6.5 * f0 + bin_range * df_amb;
        bb_hi_hz = 7.0 * f0 - bin_range * df_amb;
        bb_lo_i  = max(min(bb_lo_hz, bb_hi_hz), 1e-4) / (fs_i/2);
        bb_hi_i  = min(max(bb_lo_hz, bb_hi_hz), fs_i/2 - 1) / (fs_i/2);
        bb_hi_i  = min(bb_hi_i, 1-1e-4);
        [b_i, a_i] = butter(4, [bb_lo_i, bb_hi_i], 'bandpass');
 
        ambient_bb_std  = std(filtfilt(b_i, a_i, ambient_frame_i));
        target_bb_std_i = std(filtfilt(b_i, a_i, frame_i));
 
        % ---- Harmonics (3): H1, H2, H3 --------------------------------
        % Find index of each harmonic in f_i
        h_idx_vec = arrayfun(@(hf) find(abs(f_i - hf) == min(abs(f_i - hf)), 1), ...
                             harmonic_freqs);
        h_amp_lin = pxx_i(h_idx_vec);   % linear power at peak bin
        h_dB      = 10*log10(h_amp_lin);
 
        % ---- Ultraharmonics (3): U1, U2, U3 ---------------------------
        uh_idx_vec = zeros(1, numel(uh_orders));
        uh_amp_lin = zeros(1, numel(uh_orders));
        for k = 1:numel(uh_orders)
            f_target = uh_orders(k) * f0;
            idx_w = find(f_i >= f_target - search_range_i & ...
                         f_i <= f_target + search_range_i);
            if ~isempty(idx_w)
                [max_val, max_loc] = max(pxx_i(idx_w));
                uh_idx_vec(k) = idx_w(max_loc);
                uh_amp_lin(k) = max_val;
            end
        end
        uh_dB = 10*log10(uh_amp_lin);
 
        % ---- Threshold flags -------------------------------------------
        exceeds_harm = h_dB  > harmonic_threshold_dB;   % 1×3 logical
        exceeds_uh   = uh_dB > uh_threshold_dB;         % 1×3 logical
 
        % ---- Safety flags ----------------------------------------------
        safety_1s = target_bb_std_i > 2 * ambient_bb_std;
        safety_2s = target_bb_std_i > 3 * ambient_bb_std;
        safety_3s = target_bb_std_i > 4 * ambient_bb_std;
 
        % ---- Console output --------------------------------------------
        fprintf('  Frame %d/%d | Pressure %.4f MPa\n', f, N_frames, pressure_MPa);
        fprintf('    H_dB:  [%.2f  %.2f  %.2f]\n',  h_dB(1),  h_dB(2),  h_dB(3));
        fprintf('    UH_dB: [%.2f  %.2f  %.2f]\n', uh_dB(1), uh_dB(2), uh_dB(3));
        fprintf('    Exceeds harm: [%d %d %d]  uh: [%d %d %d]\n', ...
            exceeds_harm(1), exceeds_harm(2), exceeds_harm(3), ...
            exceeds_uh(1),   exceeds_uh(2),   exceeds_uh(3));
        fprintf('    BB std: %.4f  Ambient: %.4f  Safety: 1s=%d 2s=%d 3s=%d\n', ...
            target_bb_std_i, ambient_bb_std, safety_1s, safety_2s, safety_3s);
 
        % ---- Store as flat struct (one scalar/string per field) --------
        % Pre-define ALL fields with defaults so struct shape is identical
        % every iteration — prevents vertcat field-mismatch errors.
        R = struct( ...
            'concentration',         '', ...
            'filename',              '', ...
            'frame',                 NaN, ...
            'pressure_MPa',          NaN, ...
            'H1_idx',                NaN, ...
            'H2_idx',                NaN, ...
            'H3_idx',                NaN, ...
            'H1_dB',                 NaN, ...
            'H2_dB',                 NaN, ...
            'H3_dB',                 NaN, ...
            'U1_idx',                NaN, ...
            'U2_idx',                NaN, ...
            'U3_idx',                NaN, ...
            'U1_dB',                 NaN, ...
            'U2_dB',                 NaN, ...
            'U3_dB',                 NaN, ...
            'target_bb_std',         NaN, ...
            'ambient_bb_std',        NaN, ...
            'harmonic_threshold_dB', NaN, ...
            'U1_threshold_dB',       NaN, ...
            'U2_threshold_dB',       NaN, ...
            'U3_threshold_dB',       NaN, ...
            'exceeds_H1',            false, ...
            'exceeds_H2',            false, ...
            'exceeds_H3',            false, ...
            'exceeds_U1',            false, ...
            'exceeds_U2',            false, ...
            'exceeds_U3',            false, ...
            'safety_1s',             false, ...
            'safety_2s',             false, ...
            'safety_3s',             false  ...
        );
 
        R.concentration     = concentration;
        R.filename          = raw_name;
        R.frame             = f;
        R.pressure_MPa      = pressure_MPa;
 
        % Harmonic indices
        R.H1_idx            = h_idx_vec(1);
        R.H2_idx            = h_idx_vec(2);
        R.H3_idx            = h_idx_vec(3);
 
        % Harmonic dB
        R.H1_dB             = h_dB(1);
        R.H2_dB             = h_dB(2);
        R.H3_dB             = h_dB(3);
 
        % Ultraharmonic indices
        R.U1_idx            = uh_idx_vec(1);
        R.U2_idx            = uh_idx_vec(2);
        R.U3_idx            = uh_idx_vec(3);
 
        % Ultraharmonic dB
        R.U1_dB             = uh_dB(1);
        R.U2_dB             = uh_dB(2);
        R.U3_dB             = uh_dB(3);
 
        % Broadband
        R.target_bb_std     = target_bb_std_i;
        R.ambient_bb_std    = ambient_bb_std;
 
        % Threshold values used
        R.harmonic_threshold_dB = harmonic_threshold_dB;   % scalar
        R.U1_threshold_dB   = uh_threshold_dB(1);
        R.U2_threshold_dB   = uh_threshold_dB(2);
        R.U3_threshold_dB   = uh_threshold_dB(3);
 
        % Threshold flags
        R.exceeds_H1        = exceeds_harm(1);
        R.exceeds_H2        = exceeds_harm(2);
        R.exceeds_H3        = exceeds_harm(3);
        R.exceeds_U1        = exceeds_uh(1);
        R.exceeds_U2        = exceeds_uh(2);
        R.exceeds_U3        = exceeds_uh(3);
 
        % Safety flags
        R.safety_1s         = safety_1s;
        R.safety_2s         = safety_2s;
        R.safety_3s         = safety_3s;
 
        Results = [Results; R];
    end
end
 
%% --- Save probe CSV ---
T_probe = struct2table(Results);
writetable(T_probe, 'Probe_Features.csv');
fprintf('\nProbe CSV saved: %d rows.\n', height(T_probe));
 

%% Functions

function out = extract_probe_data(probe_path, i)
    data_files = dir(fullfile(probe_path, '*.pacq*'));
    file_path  = string(fullfile(data_files(i).folder, data_files(i).name));
    acq        = imaging.common.acquisition_reader(file_path);
    disp(acq.data_params)
    disp(acq)

    N_t      = acq.rfdata_dimensions(1);
    N_ch     = acq.rfdata_dimensions(2);
    N_frames = acq.number_of_frames;
    fs       = acq.data_params.fs;

    rf_data  = zeros(N_t, N_ch, N_frames, 'double');
    for k = 1:N_frames
        rf_data(:, :, k) = double(acq.data{k});
    end

    out.rf_data  = rf_data;
    out.t        = (0:N_t-1).' / fs;
    out.frame_t  = acq.elapsed_timevector(:);
    out.fs       = fs;
end

function features = get_features(f, pxx, f0)
    pxx_db = 10*log10(pxx);
    df     = f(2) - f(1);

    harmonics      = 6:12;
    ultraharmonics = 6.5:1:12.5;
    search_window  = round(0.08e6 / df);

    features        = struct();
    harm_freq_vals  = zeros(1, length(harmonics));
    ultra_freq_vals = zeros(1, length(ultraharmonics));
end 
