% in-vitro PCD Dataset (simplified)
clear; clc; close all;

%% === PATH SETUP ===
data_path    = '/Users/shwethasuresh/Desktop/capstone/data/1_and_2_2026-01-16_single_tube_sweep/pcd/100x_01';
ambient_path = '/Users/shwethasuresh/Desktop/capstone/data/1_and_2_2026-01-16_single_tube_sweep/pcd/ambient';

threshold_N = 1;   % threshold = mean + N * std  (across ambient files)

fprintf('Using data folder:    %s\n', data_path);
fprintf('Using ambient folder: %s\n', ambient_path);

if ~exist(data_path, 'dir')
    error('Data folder not found: %s', data_path);
end
if ~exist(ambient_path, 'dir')
    error('Ambient folder not found: %s', ambient_path);
end

%% === LOAD TARGET FILE ===
mat_files = dir(fullfile(data_path, '*.mat'));
if isempty(mat_files)
    error('No .mat files found in: %s', data_path);
end

mat_files_alt = mat_files(2:2:end);

fprintf('\nFound %d alternating .mat files:\n', length(mat_files_alt));
for i = 1:length(mat_files_alt)
    fprintf('%d. %s\n', i, mat_files_alt(i).name);
end

file_number = 32;
if file_number > length(mat_files_alt)
    error('file_number exceeds available files.');
end

file_to_load = fullfile(data_path, mat_files_alt(file_number).name);
reader = load(file_to_load);

fprintf('\nData loaded successfully\n');
fprintf('\n=== DATA INFORMATION ===\n');
fprintf('File name:            %s\n',   mat_files_alt(file_number).name);
fprintf('Tstart:               %d s\n', reader.Tstart);
fprintf('Tinterval:            %d s\n', reader.Tinterval);
fprintf('Extra Samples:        %.0f\n', reader.ExtraSamples);
fprintf('Requested Length:     %.2f mm\n', reader.RequestedLength);
fprintf('Length:               %.2f mm\n', reader.Length);
fprintf('Number of versions:   %.0f\n', length(reader.Version));
fprintf('A length:             %.0f\n', length(reader.A));
fprintf('B length:             %.0f\n', length(reader.B));

%% === TIME VECTOR ===
t = reader.Tstart + (0:length(reader.B)-1) * reader.Tinterval;

%% === WELCH PSD (target file) ===
[pxx, f] = pwelch(reader.B, [], [], [], 1/reader.Tinterval);
WelchSpec = sqrt(pxx);
df = f(2) - f(1);

%% === SIGNAL STATS ===
fprintf('\n=== SIGNAL ANALYSIS ===\n');
fprintf('  Mean:          %.4f mV\n', mean(reader.B));
fprintf('  Std deviation: %.4f mV\n', std(reader.B));
fprintf('  Peak-to-peak:  %.4f mV\n', max(reader.B) - min(reader.B));
fprintf('  RMS:           %.4f mV\n', rms(reader.B));

%% === DRIVE / HARMONIC SETTINGS ===
drive_frequency = 0.5e6;
n_harmonics     = 4;
harmonic_freqs  = drive_frequency * (1:n_harmonics);
bin_range       = 1000;

uh_orders           = [1.5 2.5 3.5];
ultraharmonic_freqs = uh_orders * drive_frequency;
search_range        = bin_range * df;

%% === HARMONICS (target file) ===
harmonic_indices = zeros(1, n_harmonics);
for k = 1:n_harmonics
    [~, idx] = min(abs(f - harmonic_freqs(k)));
    harmonic_indices(k) = idx;
end
harmonic_amplitudes = WelchSpec(harmonic_indices);

%% === ULTRAHARMONICS (target file) ===
ultraharmonic_indices    = zeros(1, length(uh_orders));
ultraharmonic_amplitudes = zeros(1, length(uh_orders));

for k = 1:length(uh_orders)
    f_target   = ultraharmonic_freqs(k);
    idx_window = find(f >= f_target - search_range & f <= f_target + search_range);
    if isempty(idx_window)
        warning('No points found in search range for U%.1f', uh_orders(k));
        continue
    end
    [max_val, max_idx]          = max(WelchSpec(idx_window));
    ultraharmonic_indices(k)    = idx_window(max_idx);
    ultraharmonic_amplitudes(k) = max_val;
end



%% === AMBIENT THRESHOLD COMPUTATION ===
% Load the SINGLE ambient file that matches file_number (same pressure level)

amb_files = dir(fullfile(ambient_path, '*.mat'));
if isempty(amb_files)
    error('No .mat files found in ambient folder: %s', ambient_path);
end

amb_files_alt = amb_files(2:2:end);
n_amb = length(amb_files_alt);

% --- Guard: ensure the matched index exists ---
if file_number > n_amb
    error(['file_number (%d) exceeds the number of alternating ambient files (%d).\n' ...
           'Cannot form a matched ambient baseline.'], file_number, n_amb);
end

matched_amb_file = fullfile(ambient_path, amb_files_alt(file_number).name);
fprintf('\nMatched ambient file (same pressure index %d): %s\n', ...
    file_number, amb_files_alt(file_number).name);

% Compute Welch spectrum for the matched ambient file
amb_data = load(matched_amb_file);
[amb_pxx, amb_f] = pwelch(amb_data.B, [], [], [], 1/amb_data.Tinterval);
amb_spec = sqrt(amb_pxx);
amb_df   = amb_f(2) - amb_f(1);
amb_sr   = bin_range * amb_df;

% UH peak amplitudes (dB) from the matched ambient file
amb_uh_peak_dB = NaN(1, length(uh_orders));
for k = 1:length(uh_orders)
    f_target   = uh_orders(k) * drive_frequency;
    idx_window = find(amb_f >= f_target - amb_sr & amb_f <= f_target + amb_sr);
    if isempty(idx_window)
        warning('No ambient points found for U%.1f', uh_orders(k));
        continue
    end
    peak_amp           = max(amb_spec(idx_window));
    amb_uh_peak_dB(k)  = 20*log10(peak_amp + eps);
end

% -----------------------------------------------------------------
% Threshold = ambient peak + N * (per-UH noise std of ambient file)
%
% "std" here is estimated from the spectral variance inside the same
% search window of the AMBIENT file, giving a per-UH spread estimate
% even from a single file.
% -----------------------------------------------------------------
amb_uh_std_dB = NaN(1, length(uh_orders));
for k = 1:length(uh_orders)
    f_target   = uh_orders(k) * drive_frequency;
    idx_window = find(amb_f >= f_target - amb_sr & amb_f <= f_target + amb_sr);
    if isempty(idx_window), continue; end
    window_vals_dB      = 20*log10(amb_spec(idx_window) + eps);
    amb_uh_std_dB(k)    = std(window_vals_dB);
end

threshold_dB = amb_uh_peak_dB + threshold_N * amb_uh_std_dB;

% Target file UH peak dB (for comparison)
target_uh_dB = 20*log10(ultraharmonic_amplitudes + eps);

% Exceeds threshold flag (UH)
exceeds = target_uh_dB > threshold_dB;

%% === HARMONIC THRESHOLD COMPUTATION ===
% For harmonics, the ambient file at the same pressure ALREADY contains
% harmonic energy from the transducer drive, so using the ambient peak
% inside the harmonic window as a baseline is too conservative.
% Instead, threshold = ambient noise floor (mean PSD in inter-harmonic
% gaps) + N * std of that noise floor — i.e. "is this harmonic above
% the ambient noise?", not "is it above ambient harmonic energy?".

% Estimate ambient noise floor from inter-harmonic regions
% Use gaps between harmonics, excluding ±bin_range bins around each harmonic
amb_noise_mask = true(size(amb_f));
for k = 1:n_harmonics
    [~, h_idx] = min(abs(amb_f - harmonic_freqs(k)));
    lo = max(1, h_idx - bin_range);
    hi = min(length(amb_f), h_idx + bin_range);
    amb_noise_mask(lo:hi) = false;
end
% Also exclude UH regions
for k = 1:length(uh_orders)
    f_target = uh_orders(k) * drive_frequency;
    [~, u_idx] = min(abs(amb_f - f_target));
    lo = max(1, u_idx - bin_range);
    hi = min(length(amb_f), u_idx + bin_range);
    amb_noise_mask(lo:hi) = false;
end

amb_noise_vals_dB  = 20*log10(amb_spec(amb_noise_mask) + eps);
amb_noise_mean_dB  = mean(amb_noise_vals_dB);
amb_noise_std_dB   = std(amb_noise_vals_dB);
harmonic_threshold_dB_scalar = amb_noise_mean_dB + threshold_N * amb_noise_std_dB;

% Same threshold applied to all harmonics (H2, H3, H4)
harmonic_threshold_dB = NaN(1, n_harmonics);
harmonic_threshold_dB(2:end) = harmonic_threshold_dB_scalar;

% For display: store ambient noise stats per harmonic column (same value repeated)
amb_harm_peak_dB = repmat(amb_noise_mean_dB, 1, n_harmonics);
amb_harm_std_dB  = repmat(amb_noise_std_dB,  1, n_harmonics);

target_harm_dB  = 20*log10(harmonic_amplitudes + eps);
% harmonic_amplitudes(1) = f0, (2) = H1, (3) = H2, (4) = H3
% harmonic_threshold_dB(1) = NaN (f0 skipped), (2:end) = computed
exceeds_harm    = false(1, n_harmonics);
exceeds_harm(1) = true;   % f0 always shown
for k = 2:n_harmonics
    exceeds_harm(k) = target_harm_dB(k) > harmonic_threshold_dB(k);
end

%% === PRINT THRESHOLD TABLES ===

% --- Harmonic table (H2, H3, H4) ---
fprintf('\n=== Harmonic Threshold Summary ===\n');
fprintf('Matched ambient file index: %d  →  %s\n', file_number, amb_files_alt(file_number).name);
fprintf('Threshold = Ambient Noise Floor Mean + %d × Std  (inter-harmonic gaps)  [f0 always shown]\n', threshold_N);
fprintf('Ambient noise floor: Mean = %.2f dB,  Std = %.2f dB,  Threshold = %.2f dB\n\n', ...
    amb_noise_mean_dB, amb_noise_std_dB, harmonic_threshold_dB_scalar);
fprintf('%-6s  %12s  %10s  %8s\n', 'H', 'Target Peak', 'Threshold', 'Exceeds?');
fprintf('%s\n', repmat('-', 1, 45));
for k = 2:n_harmonics
    flag = 'NO';
    if exceeds_harm(k), flag = 'YES ***'; end
    fprintf('H%-5d  %+12.2f  %+10.2f  %s\n', ...
        k-1, target_harm_dB(k), harmonic_threshold_dB(k), flag);
end
fprintf('%s\n', repmat('-', 1, 45));

% --- Ultraharmonic table ---
fprintf('\n=== Ultraharmonic Threshold Summary ===\n');
fprintf('Matched ambient file index: %d  →  %s\n', file_number, amb_files_alt(file_number).name);
fprintf('Threshold = Ambient Peak + %d × Std(ambient search window)\n\n', threshold_N);
fprintf('%-6s  %12s  %12s  %10s  %12s  %8s\n', ...
    'UH', 'Amb Peak (dB)', 'Amb Std (dB)', 'Threshold', 'Target Peak', 'Exceeds?');
fprintf('%s\n', repmat('-', 1, 70));
for k = 1:length(uh_orders)
    flag = 'NO';
    if exceeds(k), flag = 'YES ***'; end
    fprintf('U%-5.1f  %+12.2f  %+12.2f  %+10.2f  %+12.2f  %s\n', ...
        uh_orders(k), ...
        amb_uh_peak_dB(k), ...
        amb_uh_std_dB(k), ...
        threshold_dB(k), ...
        target_uh_dB(k), ...
        flag);
end
fprintf('%s\n', repmat('-', 1, 70));
fprintf('All values in dB (20·log10 of Welch amplitude spectrum)\n');

%% === SAFETY ASSESSMENT (U1.5 Broadband Region) ===
% Bandpass filter both target and matched ambient to the U1.5 BB region
% [1.5*f0 + bin_range*df  →  2.0*f0 - bin_range*df]

fs          = 1 / reader.Tinterval;
bb_lo       = (1.5 * drive_frequency + bin_range * df) / (fs/2);
bb_hi       = (2.0 * drive_frequency - bin_range * df) / (fs/2);

% Clamp to valid range
bb_lo = max(bb_lo, 1e-4);
bb_hi = min(bb_hi, 1 - 1e-4);

[b_bp, a_bp] = butter(4, [bb_lo, bb_hi], 'bandpass');

target_filtered  = filtfilt(b_bp, a_bp, double(reader.B));
ambient_filtered = filtfilt(b_bp, a_bp, double(amb_data.B));

target_bb_std  = std(target_filtered);
ambient_bb_std = std(ambient_filtered);
ambient_bb_mean = mean(ambient_filtered);  % should be ~0 for bandpassed signal

% Thresholds: ambient std + 1,2,3 x ambient std (i.e. 2,3,4 x amb std total)
thresh_1s = ambient_bb_std * 2;   % mean≈0, so mean+1σ = 1σ+1σ
thresh_2s = ambient_bb_std * 3;
thresh_3s = ambient_bb_std * 4;

% Determine safety level
if target_bb_std > thresh_3s
    safety_label = 'UNSAFE (> 3σ)';
elseif target_bb_std > thresh_2s
    safety_label = 'UNSAFE (> 2σ)';
elseif target_bb_std > thresh_1s
    safety_label = 'UNSAFE (> 1σ)';
else
    safety_label = 'SAFE';
end

fprintf('\n=== Safety Assessment: U1.5 Broadband Region ===\n');
fprintf('Bandpass range:        %.4f – %.4f MHz\n', ...
    (1.5*drive_frequency + bin_range*df)/1e6, ...
    (2.0*drive_frequency - bin_range*df)/1e6);
fprintf('Ambient BB std:        %.6f\n', ambient_bb_std);
fprintf('Target  BB std:        %.6f\n', target_bb_std);
fprintf('\n  Threshold 1σ:  %.6f  →  %s\n', thresh_1s, string_exceed(target_bb_std > thresh_1s));
fprintf('  Threshold 2σ:  %.6f  →  %s\n', thresh_2s, string_exceed(target_bb_std > thresh_2s));
fprintf('  Threshold 3σ:  %.6f  →  %s\n', thresh_3s, string_exceed(target_bb_std > thresh_3s));
fprintf('\n  *** OVERALL: %s ***\n', safety_label);
fprintf('%s\n', repmat('-', 1, 50));

%% === LEGACY CONSOLE OUTPUT ===
fprintf('\n=== Ultraharmonic Peak Detection ===\n');
for k = 1:length(uh_orders)
    fprintf('U%.1f: %.6f MHz, Amplitude = %.2f dB\n', ...
        uh_orders(k), f(ultraharmonic_indices(k))/1e6, ...
        20*log10(ultraharmonic_amplitudes(k)));
end

fprintf('\n=== Total Harmonic Power (dB) ===\n');
for k = 1:n_harmonics
    if k == 1, label = 'f0'; else, label = sprintf('H%d', k-1); end
    startBin = max(1, harmonic_indices(k) - bin_range);
    endBin   = min(length(pxx), harmonic_indices(k) + bin_range);
    power_dB = 10*log10(sum(pxx(startBin:endBin)) * df);
    fprintf('%s: %.3f MHz, Total Power = %.2f dB\n', ...
        label, f(harmonic_indices(k))/1e6, power_dB);
end

fprintf('\n=== Ultraharmonic Total Power (dB) ===\n');
for k = 1:length(uh_orders)
    startBin = max(1, ultraharmonic_indices(k) - bin_range);
    endBin   = min(length(pxx), ultraharmonic_indices(k) + bin_range);
    power_dB = 10*log10(sum(pxx(startBin:endBin)) * df);
    fprintf('U%.1f: %.3f MHz, Total Power = %.2f dB\n', ...
        uh_orders(k), f(ultraharmonic_indices(k))/1e6, power_dB);
end



%% === PLOT ===
f_norm      = f / drive_frequency;
spectrum_dB = 10*log10(pxx + eps);
dark_green  = [0 0.5 0];

figure('Color', 'w', 'Position', [200 200 900 600]);
ax = axes;
hold(ax, 'on');

% Full spectrum (gray)
plot(ax, f_norm, spectrum_dB, ...
    'Color', [0.75 0.75 0.75], 'LineWidth', 1.2, ...
    'DisplayName', 'Full Acoustic Spectrum');

% Broadband noise highlight — shade between each UH and the next harmonic
% bb_regions: [f_start, f_end] in normalised units (f/f0)
bin_range_norm = bin_range * df / drive_frequency;  % bin_range in normalised units
bb_regions = [1.5 + bin_range_norm, 2.0 - bin_range_norm;
              2.5 + bin_range_norm, 3.0 - bin_range_norm;
              3.5 + bin_range_norm, 4.0 - bin_range_norm];

first_plotted = false;
for k = 1:size(bb_regions, 1)
    mask = f_norm >= bb_regions(k,1) & f_norm <= bb_regions(k,2);
    if ~any(mask), continue; end
    seg = NaN(size(spectrum_dB));
    seg(mask) = spectrum_dB(mask);
    if ~first_plotted
        plot(ax, f_norm, seg, 'b', 'LineWidth', 2, ...
            'DisplayName', 'Broadband Noise Regions');
        first_plotted = true;
    else
        plot(ax, f_norm, seg, 'b', 'LineWidth', 2, ...
            'HandleVisibility', 'off');
    end
end

% Harmonics (red circles) — f0 always; H2+ only if exceeds harmonic threshold
legend_added_harm = false;
for k = 1:n_harmonics
    if ~exceeds_harm(k), continue; end
    if ~legend_added_harm
        plot(ax, f_norm(harmonic_indices(k)), ...
            20*log10(harmonic_amplitudes(k) + eps), ...
            'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, ...
            'DisplayName', 'Fundamental (f_0) and Harmonics (exceeds threshold)');
        legend_added_harm = true;
    else
        plot(ax, f_norm(harmonic_indices(k)), ...
            20*log10(harmonic_amplitudes(k) + eps), ...
            'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, ...
            'HandleVisibility', 'off');
    end
end

% Harmonic threshold lines (H2, H3, H4) — same style as UH thresholds
harm_thr_colors = [0.85 0.3 0;   % burnt orange for H2
                   0.75 0.1 0;   % darker orange for H3
                   0.60 0.0 0];  % deep red-brown for H4
harm_labels = {'H1','H2','H3'};

for k = 2:n_harmonics
    if isnan(harmonic_threshold_dB(k)), continue; end
    f_ctr  = harmonic_freqs(k) / drive_frequency;
    half_w = search_range / drive_frequency;
    x_seg  = [f_ctr - half_w, f_ctr + half_w];
    y_seg  = [harmonic_threshold_dB(k), harmonic_threshold_dB(k)];

    plot(ax, x_seg, y_seg, '--', ...
        'Color',     harm_thr_colors(k-1,:), ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('Threshold %s = %.1f dB', harm_labels{k-1}, harmonic_threshold_dB(k)));

    text(ax, f_ctr, harmonic_threshold_dB(k) + 2, ...
        sprintf('Thr %s', harm_labels{k-1}), ...
        'Color', harm_thr_colors(k-1,:), 'FontSize', 9, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

% Ultraharmonics (green triangles) — only plot if target peak exceeds threshold
legend_added_uh = false;
for k = 1:length(uh_orders)
    if ultraharmonic_indices(k) <= 0, continue; end
    if ~exceeds(k), continue; end

    if ~legend_added_uh
        plot(ax, f_norm(ultraharmonic_indices(k)), ...
            20*log10(ultraharmonic_amplitudes(k) + eps), ...
            '^', 'Color', dark_green, 'MarkerFaceColor', dark_green, 'MarkerSize', 8, ...
            'DisplayName', 'Ultraharmonics (exceeds threshold)');
        legend_added_uh = true;
    else
        plot(ax, f_norm(ultraharmonic_indices(k)), ...
            20*log10(ultraharmonic_amplitudes(k) + eps), ...
            '^', 'Color', dark_green, 'MarkerFaceColor', dark_green, 'MarkerSize', 8, ...
            'HandleVisibility', 'off');
    end
end

% Threshold lines per UH
uh_colors = [0.85 0 0.85;
             0.60 0 0.80;
             0.30 0 0.70];

for k = 1:length(uh_orders)
    f_ctr    = uh_orders(k);
    half_w   = search_range / drive_frequency;
    x_seg    = [f_ctr - half_w, f_ctr + half_w];
    y_seg    = [threshold_dB(k), threshold_dB(k)];

    plot(ax, x_seg, y_seg, '--', ...
        'Color',     uh_colors(k,:), ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('Threshold U%.1f = %.1f dB', uh_orders(k), threshold_dB(k)));

    text(ax, f_ctr, threshold_dB(k) + 2, ...
        sprintf('Thr U%.1f', uh_orders(k)), ...
        'Color', uh_colors(k,:), 'FontSize', 9, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

% Harmonic labels — f0 always; H2+ only if exceeds threshold
label_offset = 5;
for k = 1:n_harmonics
    if ~exceeds_harm(k), continue; end
    if k == 1, label = 'f_0'; else, label = sprintf('H%d', k-1); end
    text(ax, f_norm(harmonic_indices(k)), ...
        20*log10(harmonic_amplitudes(k)+eps) + label_offset, label, ...
        'Color', 'r', 'FontWeight', 'bold', 'FontSize', 12, ...
        'HorizontalAlignment', 'center');
end

% Ultraharmonic labels — only if exceeds threshold
for k = 1:length(uh_orders)
    if ultraharmonic_indices(k) <= 0, continue; end
    if ~exceeds(k), continue; end
    text(ax, f_norm(ultraharmonic_indices(k)), ...
        20*log10(ultraharmonic_amplitudes(k)+eps) + label_offset, ...
        sprintf('U%.1f', uh_orders(k)), ...
        'Color', dark_green, 'FontWeight', 'bold', 'FontSize', 12, ...
        'HorizontalAlignment', 'center');
end

% Formatting
title(ax, 'Normalized Cavitation Spectrum — Harmonics, Ultraharmonics, Broadband Noise, and Thresholds', ...
    'FontSize', 14, 'FontWeight', 'bold');
xlabel(ax, 'Normalized Frequency  (f / f_0)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel(ax, 'Acoustic Power Spectral Density (dB)', 'FontSize', 13, 'FontWeight', 'bold');
set(ax, 'FontSize', 12, 'LineWidth', 1.2);
xlim(ax, [0 6]);
ylim(ax, [-150 -50]);
grid(ax, 'on');
box(ax, 'on');
legend(ax, 'Location', 'northeast');

%% === HELPER ===
function s = string_exceed(flag)
    if flag, s = 'EXCEEDS'; else, s = 'within'; end
end


