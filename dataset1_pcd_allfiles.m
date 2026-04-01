%% in-vitro PCD Dataset — Batch Harmonic Analysis
clear; clc; close all;

%% === PATH SETUP ===
scriptDir = fileparts(mfilename('fullpath'));
cd(scriptDir);
fprintf('Script directory: %s\n', scriptDir);

repo_path = fullfile(scriptDir, '04_Projects', '17_FUS_Instruments_Interface');
if exist(repo_path, 'dir')
    addpath(genpath(repo_path));
    fprintf('Added repository to path\n');
else
    warning('Repo not found — continuing anyway');
end

pcd_base     = '/Users/shwethasuresh/Desktop/capstone/data/1_and_2_2026-01-16_single_tube_sweep/pcd';
ambient_path = fullfile(pcd_base, 'ambient');
assert(exist(ambient_path, 'dir') == 7, 'Ambient folder not found');

data_folders = { ...
    fullfile(pcd_base, '1m_x_01');   ...
    fullfile(pcd_base, '1m_x_02');   ...
    fullfile(pcd_base, '100x_01');   ...
    fullfile(pcd_base, '100x_02');   ...
    fullfile(pcd_base, '1000x_01');  ...
    fullfile(pcd_base, '1000x_02') };

% Verify all folders exist
for d = 1:length(data_folders)
    assert(exist(data_folders{d}, 'dir') == 7, 'Data folder not found: %s', data_folders{d});
end

%% === ANALYSIS PARAMETERS ===
params.drive_frequency  = 0.5e6;
params.n_harmonics      = 4;
params.uh_orders        = [1.5 2.5 3.5];
params.bin_range        = 1000;
params.noise_exclusion  = 50;
params.uh_upper_harm    = [2, 3, 4];
params.threshold_N      = 1;   % mean + N*std

%% === LOAD FILE LISTS ===
amb_files     = dir(fullfile(ambient_path, '*.mat'));
amb_files_alt = amb_files(2:2:end);
N_amb         = length(amb_files_alt);
assert(N_amb > 0, 'No ambient .mat files found');

% Collect all target files across all folders, preserving folder label
% and per-folder index (used for matched ambient lookup)
all_target_files  = struct('name', {}, 'folder', {}, 'folder_label', {}, 'folder_index', {});
for d = 1:length(data_folders)
    mf  = dir(fullfile(data_folders{d}, '*.mat'));
    mfa = mf(2:2:end);
    [~, folder_label] = fileparts(data_folders{d});
    for j = 1:length(mfa)
        entry.name         = mfa(j).name;
        entry.folder       = mfa(j).folder;
        entry.folder_label = folder_label;
        entry.folder_index = j;   % per-folder index → matches ambient file j
        all_target_files(end+1) = entry; 
    end
end

N = length(all_target_files);
fprintf('Found %d total target files across %d folders, %d ambient files\n', ...
    N, length(data_folders), N_amb);
assert(N > 0, 'No target .mat files found');

% No global N truncation needed — each file uses its own folder_index
% which is checked individually against N_amb inside the loop

%% === PREALLOCATE RESULTS ===
Results(N) = struct( ...
    'filename',               [], ...
    'folder',                 [], ...
    'folder_label',           [], ...
    'harmonic_power',         [], ...
    'harmonic_power_dB',      [], ...
    'ultraharmonic_power',    [], ...
    'ultraharmonic_power_dB', [], ...
    'harmonic_indices',       [], ...
    'ultraharmonic_indices',  [], ...
    'broadband_noise_power',    [], ...
    'broadband_noise_power_dB', [], ...
    'noise_percentage',         [], ...
    'exceeds_harm',           [], ...
    'exceeds_uh',             [], ...
    'harmonic_threshold_dB',  [], ...
    'uh_threshold_dB',        [], ...
    'safety_1s',              [], ...
    'safety_2s',              [], ...
    'safety_3s',              [], ...
    'target_bb_std',          [], ...
    'ambient_bb_std',         [] );

%% === MAIN LOOP ===
for i = 1:N
    fprintf('\n[%d/%d] [%s] %s\n', i, N, all_target_files(i).folder_label, all_target_files(i).name);

    % --- Load target ---
    reader = load(fullfile(all_target_files(i).folder, all_target_files(i).name));
    [pxx, f] = compute_pwelch(reader.B, reader.Tinterval);
    df = f(2) - f(1);

    amb_idx = all_target_files(i).folder_index;
    if amb_idx > N_amb
        warning('[%d/%d] folder_index %d exceeds ambient file count (%d), skipping.', i, N, amb_idx, N_amb);
        continue
    end

    % --- Load matched ambient (same pressure index within this folder) ---
    amb_data = load(fullfile(ambient_path, amb_files_alt(amb_idx).name));
    [amb_pxx, amb_f] = compute_pwelch(amb_data.B, amb_data.Tinterval);
    amb_spec = sqrt(amb_pxx);
    amb_df   = amb_f(2) - amb_f(1);
    amb_sr   = params.bin_range * amb_df;

    % --- Find harmonics & ultraharmonics ---
    [harm_idx, uh_idx] = find_harmonics(f, pxx, params);

    % --- Integrate power ---
    [Hpow, Hpow_dB] = integrate_band_power(pxx, f, harm_idx, params.bin_range);
    [Upow, Upow_dB] = integrate_band_power(pxx, f, uh_idx,   params.bin_range);

    % --- Broadband noise ---
    [BBnoise, BBnoise_dB, noise_pct] = calculate_broadband_noise(pxx, f, harm_idx, uh_idx, params);

    % --- UH threshold: ambient peak + N*std in ambient search window ---
    amb_uh_peak_dB = NaN(1, length(params.uh_orders));
    amb_uh_std_dB  = NaN(1, length(params.uh_orders));
    for k = 1:length(params.uh_orders)
        f_target   = params.uh_orders(k) * params.drive_frequency;
        idx_win    = find(amb_f >= f_target - amb_sr & amb_f <= f_target + amb_sr);
        if isempty(idx_win), continue; end
        vals_dB           = 20*log10(amb_spec(idx_win) + eps);
        amb_uh_peak_dB(k) = max(vals_dB);
        amb_uh_std_dB(k)  = std(vals_dB);
    end
    uh_threshold_dB = amb_uh_peak_dB + params.threshold_N * amb_uh_std_dB;
    target_uh_dB    = 20*log10(sqrt(pxx(uh_idx)) + eps);
    exceeds_uh      = target_uh_dB > uh_threshold_dB;

    % --- Harmonic threshold: ambient noise floor (inter-harmonic gaps) ---
    harmonic_freqs = params.drive_frequency * (1:params.n_harmonics);
    amb_noise_mask = true(size(amb_f));
    for k = 1:params.n_harmonics
        [~, h_idx] = min(abs(amb_f - harmonic_freqs(k)));
        lo = max(1, h_idx - params.bin_range);
        hi = min(length(amb_f), h_idx + params.bin_range);
        amb_noise_mask(lo:hi) = false;
    end
    for k = 1:length(params.uh_orders)
        [~, u_idx] = min(abs(amb_f - params.uh_orders(k)*params.drive_frequency));
        lo = max(1, u_idx - params.bin_range);
        hi = min(length(amb_f), u_idx + params.bin_range);
        amb_noise_mask(lo:hi) = false;
    end
    amb_noise_dB        = 20*log10(amb_spec(amb_noise_mask) + eps);
    harm_thresh_scalar  = mean(amb_noise_dB) + params.threshold_N * std(amb_noise_dB);
    harm_threshold_dB   = [NaN, repmat(harm_thresh_scalar, 1, params.n_harmonics-1)];

    target_harm_dB = Hpow_dB;   % already in dB per harmonic
    exceeds_harm   = false(1, params.n_harmonics);
    exceeds_harm(1) = true;     % f0 always shown
    for k = 2:params.n_harmonics
        exceeds_harm(k) = target_harm_dB(k) > harm_threshold_dB(k);
    end

    % --- Safety assessment: bandpass std of U1.5 BB region ---
    fs    = 1 / reader.Tinterval;
    bb_lo = (1.5*params.drive_frequency + params.bin_range*df) / (fs/2);
    bb_hi = (2.0*params.drive_frequency - params.bin_range*df) / (fs/2);
    bb_lo = max(bb_lo, 1e-4);
    bb_hi = min(bb_hi, 1 - 1e-4);

    [b_bp, a_bp]     = butter(4, [bb_lo, bb_hi], 'bandpass');
    target_bb_std    = std(filtfilt(b_bp, a_bp, double(reader.B)));
    ambient_bb_std   = std(filtfilt(b_bp, a_bp, double(amb_data.B)));

    thresh_1s = ambient_bb_std * 2;
    thresh_2s = ambient_bb_std * 3;
    thresh_3s = ambient_bb_std * 4;

    safety_1s = target_bb_std > thresh_1s;
    safety_2s = target_bb_std > thresh_2s;
    safety_3s = target_bb_std > thresh_3s;

    fprintf('  Safety: 1σ=%s  2σ=%s  3σ=%s  (target std=%.4e, amb std=%.4e)\n', ...
        tf2str(safety_1s), tf2str(safety_2s), tf2str(safety_3s), ...
        target_bb_std, ambient_bb_std);

    % --- Store ---
    Results(i).filename               = all_target_files(i).name;
    Results(i).folder                 = all_target_files(i).folder;
    Results(i).folder_label           = all_target_files(i).folder_label;
    Results(i).harmonic_power         = Hpow;
    Results(i).harmonic_power_dB      = Hpow_dB;
    Results(i).ultraharmonic_power    = Upow;
    Results(i).ultraharmonic_power_dB = Upow_dB;
    Results(i).harmonic_indices       = harm_idx;
    Results(i).ultraharmonic_indices  = uh_idx;
    Results(i).broadband_noise_power    = BBnoise;
    Results(i).broadband_noise_power_dB = BBnoise_dB;
    Results(i).noise_percentage         = noise_pct;
    Results(i).exceeds_harm           = exceeds_harm;
    Results(i).exceeds_uh             = exceeds_uh;
    Results(i).harmonic_threshold_dB  = harm_threshold_dB;
    Results(i).uh_threshold_dB        = uh_threshold_dB;
    Results(i).safety_1s              = safety_1s;
    Results(i).safety_2s              = safety_2s;
    Results(i).safety_3s              = safety_3s;
    Results(i).target_bb_std          = target_bb_std;
    Results(i).ambient_bb_std         = ambient_bb_std;
end

%% === RESULTS TABLE ===
H  = vertcat(Results.harmonic_power_dB);
U  = vertcat(Results.ultraharmonic_power_dB);
Hi = vertcat(Results.harmonic_indices);
Ui = vertcat(Results.ultraharmonic_indices);
BB = vertcat(Results.broadband_noise_power_dB);

% Build exceeds columns explicitly as Nx1 vectors
EH2 = false(N,1); EH3 = false(N,1); EH4 = false(N,1);
EU1 = false(N,1); EU2 = false(N,1); EU3 = false(N,1);
for i = 1:N
    eh = Results(i).exceeds_harm;
    eu = Results(i).exceeds_uh;
    if numel(eh) >= 4, EH2(i)=eh(2); EH3(i)=eh(3); EH4(i)=eh(4); end
    if numel(eu) >= 3, EU1(i)=eu(1); EU2(i)=eu(2); EU3(i)=eu(3); end
end

ResultTable = table( ...
    {Results.folder_label}', {Results.filename}', ...
    Hi(:,1), Hi(:,2), Hi(:,3), Hi(:,4), ...
    H(:,1),  H(:,2),  H(:,3),  H(:,4), ...
    Ui(:,1), Ui(:,2), Ui(:,3), ...
    U(:,1),  U(:,2),  U(:,3), ...
    BB(:,1), BB(:,2), BB(:,3), ...
    EH2, EH3, EH4, ...
    EU1, EU2, EU3, ...
    [Results.safety_1s]', ...
    [Results.safety_2s]', ...
    [Results.safety_3s]', ...
    [Results.target_bb_std]', ...
    [Results.ambient_bb_std]', ...
    'VariableNames', { ...
        'Folder', 'File', ...
        'f0_idx','H1_idx','H2_idx','H3_idx', ...
        'f0_dB','H1_dB','H2_dB','H3_dB', ...
        'U1_idx','U2_idx','U3_idx', ...
        'U1_dB','U2_dB','U3_dB', ...
        'U1_BB_dB','U2_BB_dB','U3_BB_dB', ...
        'H1_exceeds','H2_exceeds','H3_exceeds', ...
        'U1_exceeds','U2_exceeds','U3_exceeds', ...
        'Unsafe_1s','Unsafe_2s','Unsafe_3s', ...
        'Target_BB_std', 'Ambient_BB_std'} );

disp(ResultTable);
writetable(ResultTable, 'PCD_Harmonic_Features.csv');
save('PCD_Features_AllFiles.mat', 'Results', 'ResultTable');
fprintf('\n✓ PCD batch analysis complete\n');

%% LOCAL FUNCTIONS

% -----------------------------------------------------------------
% 0 = safe, 1 = unsafe
function s = tf2str(flag)
    if flag, s = 'UNSAFE'; else, s = 'SAFE'; end
end

function [pxx, f] = compute_pwelch(signal, Tinterval)
    [pxx, f] = pwelch(signal, [], [], [], 1/Tinterval);
end

% -----------------------------------------------------------------
function [harmonic_indices, ultraharmonic_indices] = find_harmonics(f, pxx, params)
    WelchSpec = sqrt(pxx);
    df        = f(2) - f(1);

    harmonic_freqs   = params.drive_frequency * (1:params.n_harmonics);
    harmonic_indices = zeros(1, params.n_harmonics);
    for k = 1:params.n_harmonics
        [~, harmonic_indices(k)] = min(abs(f - harmonic_freqs(k)));
    end

    uh_freqs              = params.uh_orders * params.drive_frequency;
    ultraharmonic_indices = zeros(1, length(uh_freqs));
    search_range          = params.bin_range * df;
    for k = 1:length(uh_freqs)
        idx_window = find(f >= uh_freqs(k) - search_range & ...
                          f <= uh_freqs(k) + search_range);
        if isempty(idx_window)
            ultraharmonic_indices(k) = NaN; continue
        end
        [~, max_idx]             = max(WelchSpec(idx_window));
        ultraharmonic_indices(k) = idx_window(max_idx);
    end
end

% -----------------------------------------------------------------
function [power_lin, power_dB] = integrate_band_power(pxx, f, indices, bin_range)
    df        = f(2) - f(1);
    power_lin = zeros(1, length(indices));
    power_dB  = zeros(1, length(indices));
    for k = 1:length(indices)
        if isnan(indices(k))
            power_lin(k) = NaN; power_dB(k) = NaN; continue
        end
        startBin     = max(1,           indices(k) - bin_range);
        endBin       = min(length(pxx), indices(k) + bin_range);
        power_lin(k) = sum(pxx(startBin:endBin)) * df;
        power_dB(k)  = 10*log10(power_lin(k));
    end
end

% -----------------------------------------------------------------
function [noise_power, noise_power_dB, noise_percentage] = ...
        calculate_broadband_noise(pxx, f, harm_idx, uh_idx, params)
    df              = f(2) - f(1);
    nUH             = length(uh_idx);
    noise_exclusion = params.noise_exclusion;

    noise_power      = zeros(1, nUH);
    noise_power_dB   = zeros(1, nUH);
    noise_percentage = zeros(1, nUH);

    for k = 1:nUH
        uh_center = uh_idx(k);
        h_hi_idx  = harm_idx(params.uh_upper_harm(k));

        if isnan(uh_center) || uh_center <= 0
            noise_power(k) = NaN; noise_power_dB(k) = NaN; continue
        end

        bb_start = max(1,           uh_center + noise_exclusion + 1);
        bb_end   = min(length(pxx), h_hi_idx  - noise_exclusion - 1);

        if bb_end <= bb_start
            warning('Empty BB gap for U%.1f', params.uh_orders(k));
            noise_power(k) = NaN; noise_power_dB(k) = NaN; continue
        end

        noise_psd            = pxx(bb_start:bb_end);
        noise_power(k)       = mean(noise_psd * df);
        noise_power_dB(k)    = 10*log10(noise_power(k));
        noise_percentage(k)  = 100 * (bb_end - bb_start + 1) / (2*params.bin_range + 1);
    end
end
