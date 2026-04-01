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

%% === BROADBAND NOISE (target file) ===
noise_exclusion  = 1000;
uh_upper_harm    = [2, 3, 4];

broadband_noise_power    = zeros(1, length(uh_orders));
broadband_noise_power_dB = zeros(1, length(uh_orders));
bb_bounds                = zeros(length(uh_orders), 2);

for k = 1:length(uh_orders)
    uh_idx   = ultraharmonic_indices(k);
    h_hi_idx = harmonic_indices(uh_upper_harm(k));
    bb_start = max(1,           uh_idx   + noise_exclusion + 1);
    bb_end   = min(length(pxx), h_hi_idx - noise_exclusion - 1);
    if bb_end <= bb_start
        warning('Empty BB region for U%.1f', uh_orders(k));
        continue
    end
    noise_psd                    = pxx(bb_start:bb_end);
    broadband_noise_power(k)     = mean(noise_psd * df);
    broadband_noise_power_dB(k)  = 10*log10(broadband_noise_power(k));
    bb_bounds(k,:)               = [bb_start, bb_end];
end

noise_floor_dB = mean(broadband_noise_power_dB, 'omitnan');

%% === AMBIENT THRESHOLD COMPUTATION ===

amb_files = dir(fullfile(ambient_path, '*.mat'));
if isempty(amb_files)
    error('No .mat files found in ambient folder: %s', ambient_path);
end

amb_files_alt = amb_files(2:2:end);
n_amb = length(amb_files_alt);
fprintf('\nAmbient files used for threshold (%d total):\n', n_amb);
for i = 1:n_amb
    fprintf('  %d. %s\n', i, amb_files_alt(i).name);
end

% Collect per-file UH peak dB values from every ambient file
% Rows = ambient files, Cols = UH order (1.5, 2.5, 3.5)
amb_uh_dB = NaN(n_amb, length(uh_orders));

for i = 1:n_amb
    try
        amb_data = load(fullfile(ambient_path, amb_files_alt(i).name));
        [amb_pxx, amb_f] = pwelch(amb_data.B, [], [], [], 1/amb_data.Tinterval);
        amb_spec = sqrt(amb_pxx);
        amb_df   = amb_f(2) - amb_f(1);
        amb_sr   = bin_range * amb_df;

        for k = 1:length(uh_orders)
            f_target   = uh_orders(k) * drive_frequency;
            idx_window = find(amb_f >= f_target - amb_sr & amb_f <= f_target + amb_sr);
            if isempty(idx_window), continue; end
            peak_amp       = max(amb_spec(idx_window));
            amb_uh_dB(i,k) = 20*log10(peak_amp + eps);
        end
    catch ME
        warning('Could not process ambient file %s: %s', amb_files_alt(i).name, ME.message);
    end
end

% Per-UH statistics across ambient files
amb_mean_dB = mean(amb_uh_dB, 1, 'omitnan');
amb_std_dB  = std( amb_uh_dB, 0, 1, 'omitnan');
threshold_dB = amb_mean_dB + threshold_N * amb_std_dB;

% Target file UH peak dB (for comparison in table)
target_uh_dB = 20*log10(ultraharmonic_amplitudes + eps);

% Exceeds threshold flag
exceeds = target_uh_dB > threshold_dB;

%% === PRINT THRESHOLD TABLE ===
fprintf('\n=== Ultraharmonic Threshold Summary (Mean + %d\xd7Std of Ambient) ===\n', threshold_N);
fprintf('%-6s  %10s  %10s  %10s  %12s  %8s\n', ...
    'UH', 'Ambient Mean', 'Ambient Std', 'Threshold', 'Target Peak', 'Exceeds?');
fprintf('%s\n', repmat('-', 1, 65));
for k = 1:length(uh_orders)
    flag = 'NO';
    if exceeds(k), flag = 'YES ***'; end
    fprintf('U%-5.1f  %+10.2f  %+10.2f  %+10.2f  %+12.2f  %s\n', ...
        uh_orders(k), ...
        amb_mean_dB(k), ...
        amb_std_dB(k), ...
        threshold_dB(k), ...
        target_uh_dB(k), ...
        flag);
end
fprintf('%s\n', repmat('-', 1, 65));
fprintf('All values in dB (20\xb7log10 of Welch amplitude spectrum)\n');
fprintf('Threshold = Ambient Mean + %d \xd7 Ambient Std\n', threshold_N);

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

fprintf('\n=== Broadband Noise per Ultraharmonic (UH right edge to next harmonic) ===\n');
for k = 1:length(uh_orders)
    bb_start = bb_bounds(k,1);
    bb_end   = bb_bounds(k,2);
    if bb_end <= bb_start, continue; end
    fprintf('U%.1f to H%d: [%.3f - %.3f] f0, BB = %.2f dB\n', ...
        uh_orders(k), uh_upper_harm(k)-1, ...
        f(bb_start)/drive_frequency, f(bb_end)/drive_frequency, ...
        broadband_noise_power_dB(k));
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

% Broadband noise highlight
first_plotted = false;
for k = 1:length(uh_orders)
    bb_start = bb_bounds(k,1);
    bb_end   = bb_bounds(k,2);
    if bb_end <= bb_start, continue; end
    seg = NaN(size(spectrum_dB));
    seg(bb_start:bb_end) = spectrum_dB(bb_start:bb_end);
    if ~first_plotted
        plot(ax, f_norm, seg, 'b', 'LineWidth', 2, ...
            'DisplayName', 'Broadband Noise (1.5f_0 to 3f_0)');
        first_plotted = true;
    else
        plot(ax, f_norm, seg, 'b', 'LineWidth', 2, ...
            'HandleVisibility', 'off');
    end
end

% Noise floor line
yline(ax, noise_floor_dB, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Broadband Noise Floor = %.2f dB', noise_floor_dB));

% Harmonics (red circles)
plot(ax, f_norm(harmonic_indices), ...
    20*log10(harmonic_amplitudes + eps), ...
    'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, ...
    'DisplayName', 'Fundamental (f_0) and Harmonics');

% Ultraharmonics (green triangles) — only plot if target peak exceeds threshold
legend_added_uh = false;
for k = 1:length(uh_orders)
    if ultraharmonic_indices(k) <= 0, continue; end
    if ~exceeds(k), continue; end   % skip if below threshold

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

% --- Threshold lines per UH (one magenta dashed line per UH at its threshold dB) ---
% Draw each threshold line only over the UH search window so they don't
% clutter the whole plot.  A short horizontal segment centred on the UH
% target frequency with width = search_range on each side.
uh_colors = [0.85 0 0.85;   % magenta-ish for U1.5
             0.60 0 0.80;   % purple-ish  for U2.5
             0.30 0 0.70];  % deep purple for U3.5

for k = 1:length(uh_orders)
    f_ctr    = uh_orders(k);                    % in normalised units
    half_w   = search_range / drive_frequency;  % convert Hz -> normalised
    x_seg    = [f_ctr - half_w, f_ctr + half_w];
    y_seg    = [threshold_dB(k), threshold_dB(k)];

    plot(ax, x_seg, y_seg, '--', ...
        'Color',     uh_colors(k,:), ...
        'LineWidth', 1.8, ...
        'DisplayName', sprintf('Threshold U%.1f = %.1f dB', uh_orders(k), threshold_dB(k)));

    % Small annotation above the line
    text(ax, f_ctr, threshold_dB(k) + 2, ...
        sprintf('Thr U%.1f', uh_orders(k)), ...
        'Color', uh_colors(k,:), 'FontSize', 9, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

% Harmonic labels
label_offset = 5;
for k = 1:n_harmonics
    if k == 1, label = 'f_0'; else, label = sprintf('H%d', k-1); end
    text(ax, f_norm(harmonic_indices(k)), ...
        20*log10(harmonic_amplitudes(k)+eps) + label_offset, label, ...
        'Color', 'r', 'FontWeight', 'bold', 'FontSize', 12, ...
        'HorizontalAlignment', 'center');
end

% Ultraharmonic labels — only if exceeds threshold
for k = 1:length(uh_orders)
    if ultraharmonic_indices(k) <= 0, continue; end
    if ~exceeds(k), continue; end   % skip if below threshold
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


