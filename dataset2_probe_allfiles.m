%% Dataset 2 Analysis 
% Capstone II - Team E(Chidna)
% JiaJia Fu 
% 4/3/26

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

% Load ambient once — shared across all target files
ambient_data  = extract_probe_data(probe_path, ambient_file_idx);
a_fs          = ambient_data.fs;
ambient_frame = ambient_data.rf_data(:, 64, 1);   % matches single-frame
[pxx_amb, f_amb] = pwelch(ambient_frame, [], [], [], a_fs);
amb_spec      = pxx_amb;
amb_df        = f_amb(2) - f_amb(1);
amb_sr        = bin_range * amb_df;

% Precompute ambient UH threshold — same method as single-frame
amb_uh_peak_dB = NaN(1, length(uh_orders));
amb_uh_std_dB  = NaN(1, length(uh_orders));
for k = 1:length(uh_orders)
    f_target   = uh_orders(k) * f0;
    idx_window = find(f_amb >= f_target - amb_sr & f_amb <= f_target + amb_sr);
    if isempty(idx_window), continue; end
    amb_uh_peak_dB(k) = 10*log10(max(amb_spec(idx_window)));
    amb_uh_std_dB(k)  = std(10*log10(amb_spec(idx_window) + eps));
end
uh_threshold_dB = amb_uh_peak_dB + threshold_N * amb_uh_std_dB;

% Precompute ambient harmonic noise floor — same method as single-frame
harmonic_freqs = f0 * harmonics;
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
amb_noise_vals_dB            = 10*log10(amb_spec(amb_noise_mask));
amb_noise_mean_dB            = mean(amb_noise_vals_dB);
amb_noise_std_dB             = std(amb_noise_vals_dB);
harmonic_threshold_dB_scalar = amb_noise_mean_dB + threshold_N * amb_noise_std_dB;
harmonic_threshold_dB        = repmat(harmonic_threshold_dB_scalar, 1, n_harmonics);

% Precompute ambient BB std for safety — same method as single-frame
bb_lo_hz = 6.5*f0 + bin_range*amb_df;
bb_hi_hz = 7.0*f0 - bin_range*amb_df;
bb_lo_n  = max(min(bb_lo_hz, bb_hi_hz), 1e-4) / (a_fs/2);
bb_hi_n  = min(max(bb_lo_hz, bb_hi_hz), a_fs/2 - 1) / (a_fs/2);
bb_hi_n  = min(bb_hi_n, 1-1e-4);
[b_bp, a_bp]   = butter(4, [bb_lo_n, bb_hi_n], 'bandpass');
ambient_bb_std = std(filtfilt(b_bp, a_bp, double(ambient_frame)));

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

Results = [];  % start empty, will append for each frame

for i = 1:N_files
    if i == ambient_file_idx
        fprintf('Skipping ambient file: %s\n', target_files(i).name);
        continue
    end

    fprintf('\n[%d/%d] Processing: %s\n', i, N_files, target_files(i).name);

    % Load probe data
    clean_data = extract_probe_data(probe_path, i);
    fs_i       = clean_data.fs;
    N_frames   = size(clean_data.rf_data, 3);

    % Precompute BB filter once per file
    bb_lo_i = max(min(bb_lo_hz, bb_hi_hz), 1e-4) / (fs_i/2);
    bb_hi_i = min(max(bb_lo_hz, bb_hi_hz), fs_i/2 - 1) / (fs_i/2);
    bb_hi_i = min(bb_hi_i, 1-1e-4);
    [b_i, a_i] = butter(4, [bb_lo_i, bb_hi_i], 'bandpass');

    for f = 1:N_frames
        frame_i = double(clean_data.rf_data(:, 64, f));
        [pxx_i, f_i] = pwelch(frame_i, [], [], [], fs_i);
        df_i = f_i(2) - f_i(1);
        search_range_i = bin_range * df_i;

        % --- Harmonics ---
        harmonic_indices_i = zeros(1, n_harmonics);
        harmonic_amplitudes_i = zeros(1, n_harmonics);
        for k = 1:n_harmonics
            [~, idx] = min(abs(f_i - harmonic_freqs(k)));
            harmonic_indices_i(k) = idx;
            harmonic_amplitudes_i(k) = pxx_i(idx);
        end
        target_harm_dB_i = 10*log10(harmonic_amplitudes_i);

        % --- Ultraharmonics ---
        ultraharmonic_indices_i = zeros(1, length(uh_orders));
        ultraharmonic_amplitudes_i = zeros(1, length(uh_orders));
        for k = 1:length(uh_orders)
            f_target = uh_orders(k) * f0;
            idx_window = find(f_i >= f_target - search_range_i & f_i <= f_target + search_range_i);
            if isempty(idx_window), continue; end
            [max_val, max_idx] = max(pxx_i(idx_window));
            ultraharmonic_indices_i(k) = idx_window(max_idx);
            ultraharmonic_amplitudes_i(k) = max_val;
        end
        target_uh_dB_i = 10*log10(ultraharmonic_amplitudes_i);

        % --- Threshold comparisons ---
        exceeds_harm_i = target_harm_dB_i > harmonic_threshold_dB;
        exceeds_uh_i   = target_uh_dB_i > uh_threshold_dB;

        % --- BB safety ---
        target_bb_std_i = std(filtfilt(b_i, a_i, frame_i));
        safety_1s_i = target_bb_std_i > 2 * ambient_bb_std;
        safety_2s_i = target_bb_std_i > 3 * ambient_bb_std;
        safety_3s_i = target_bb_std_i > 4 * ambient_bb_std;

        % --- Append to Results ---
        R = struct();
        R.filename                 = target_files(i).name;
        R.frame                    = f;
        R.harmonic_amplitudes      = harmonic_amplitudes_i;
        R.target_harm_dB           = target_harm_dB_i;
        R.ultraharmonic_amplitudes = ultraharmonic_amplitudes_i;
        R.target_uh_dB             = target_uh_dB_i;
        R.harmonic_indices         = harmonic_indices_i;
        R.ultraharmonic_indices    = ultraharmonic_indices_i;
        R.exceeds_harm             = exceeds_harm_i;
        R.exceeds_uh               = exceeds_uh_i;
        R.harmonic_threshold_dB    = harmonic_threshold_dB;
        R.uh_threshold_dB          = uh_threshold_dB;
        R.target_bb_std            = target_bb_std_i;
        R.ambient_bb_std           = ambient_bb_std;
        R.safety_1s                = safety_1s_i;
        R.safety_2s                = safety_2s_i;
        R.safety_3s                = safety_3s_i;

        Results = [Results; R];  % append per frame
    end
end

%% --- Convert Results to CSV-safe table ---
N = numel(Results);
tableStruct = struct();
fields = fieldnames(Results);

for f = 1:numel(fields)
    field = fields{f};
    vals  = {Results.(field)};
    for k = 1:N
        if isnumeric(vals{k}) || islogical(vals{k})
            vals{k} = mat2str(double(vals{k}));
        elseif ischar(vals{k}) || isstring(vals{k})
            vals{k} = string(vals{k});
        else
            vals{k} = string(vals{k});
        end
    end
    tableStruct.(field) = vals(:);  % ensure column
end

T = struct2table(tableStruct);
writetable(T, 'ML_ready_probe_data.csv');
fprintf('\nCSV saved successfully with %d rows.\n', N);

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