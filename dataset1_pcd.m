% in-vitro PCD Dataset (simplified)
clear; clc; close all;

%% === PATH SETUP ===
data_path = '/Users/shwethasuresh/Desktop/capstone/data/1_and_2_2026-01-16_single_tube_sweep/pcd/100x_01';

fprintf('Using data folder: %s\n', data_path);

if ~exist(data_path, 'dir')
    error('Data folder not found: %s', data_path);
end

%% === LOAD FILES ===
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
fprintf('File name: %s\n',             mat_files_alt(file_number).name);
fprintf('Tstart: %d s\n',              reader.Tstart);
fprintf('Tinterval: %d s\n',           reader.Tinterval);
fprintf('Extra Samples: %.0f\n',       reader.ExtraSamples);
fprintf('Requested Length: %.2f mm\n', reader.RequestedLength);
fprintf('Length: %.2f mm\n',           reader.Length);
fprintf('Number of versions: %.0f\n',  length(reader.Version));
fprintf('A length: %.0f\n',            length(reader.A));
fprintf('B length: %.0f\n',            length(reader.B));

%% === TIME VECTOR ===
t = reader.Tstart + (0:length(reader.B)-1) * reader.Tinterval;

%% === WELCH PSD ===
[pxx, f] = pwelch(reader.B, [], [], [], 1/reader.Tinterval);
WelchSpec = sqrt(pxx);
df = f(2) - f(1);

%% === SIGNAL STATS ===
fprintf('\n=== SIGNAL ANALYSIS ===\n');
fprintf('  Mean:          %.4f mV\n', mean(reader.B));
fprintf('  Std deviation: %.4f mV\n', std(reader.B));
fprintf('  Peak-to-peak:  %.4f mV\n', max(reader.B) - min(reader.B));
fprintf('  RMS:           %.4f mV\n', rms(reader.B));

%% === HARMONICS ===
drive_frequency = 0.5e6;
n_harmonics     = 4;
harmonic_freqs  = drive_frequency * (1:n_harmonics);
bin_range       = 1000;

harmonic_indices = zeros(1, n_harmonics);
for k = 1:n_harmonics
    [~, idx] = min(abs(f - harmonic_freqs(k)));
    harmonic_indices(k) = idx;
end
harmonic_amplitudes = WelchSpec(harmonic_indices);

%% === ULTRAHARMONICS (peak search within ±bin_range bins) ===
uh_orders           = [1.5 2.5 3.5];
ultraharmonic_freqs = uh_orders * drive_frequency;
search_range        = bin_range * df;

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

fprintf('\n=== Ultraharmonic Peak Detection ===\n');
for k = 1:length(uh_orders)
    fprintf('U%.1f: %.6f MHz, Amplitude = %.2f dB\n', ...
        uh_orders(k), f(ultraharmonic_indices(k))/1e6, ...
        20*log10(ultraharmonic_amplitudes(k)));
end

%% === HARMONIC POWER ===
fprintf('\n=== Total Harmonic Power (dB) ===\n');
for k = 1:n_harmonics
    if k == 1, label = 'f0'; else, label = sprintf('H%d', k-1); end
    startBin = max(1, harmonic_indices(k) - bin_range);
    endBin   = min(length(pxx), harmonic_indices(k) + bin_range);
    power_dB = 10*log10(sum(pxx(startBin:endBin)) * df);
    fprintf('%s: %.3f MHz, Total Power = %.2f dB\n', ...
        label, f(harmonic_indices(k))/1e6, power_dB);
end

%% === ULTRAHARMONIC POWER ===
fprintf('\n=== Ultraharmonic Total Power (dB) ===\n');
for k = 1:length(uh_orders)
    startBin = max(1, ultraharmonic_indices(k) - bin_range);
    endBin   = min(length(pxx), ultraharmonic_indices(k) + bin_range);
    power_dB = 10*log10(sum(pxx(startBin:endBin)) * df);
    fprintf('U%.1f: %.3f MHz, Total Power = %.2f dB\n', ...
        uh_orders(k), f(ultraharmonic_indices(k))/1e6, power_dB);
end

%% === BROADBAND NOISE (gap from each UH right edge to next harmonic left edge) ===
% For each UH, the BB region is the gap between the UH peak and the
% harmonic immediately to its right:
%   U1.5 → H1  (harmonic_indices(2))
%   U2.5 → H2  (harmonic_indices(3))
%   U3.5 → H3  (harmonic_indices(4))

noise_exclusion = 1000;

% Upper harmonic for each UH
uh_upper_harm = [2, 3, 4];  % indices into harmonic_indices

broadband_noise_power    = zeros(1, length(uh_orders));
broadband_noise_power_dB = zeros(1, length(uh_orders));
bb_bounds                = zeros(length(uh_orders), 2);  % [start, end] per UH

fprintf('\n=== Broadband Noise per Ultraharmonic (UH right edge → next harmonic) ===\n');
for k = 1:length(uh_orders)
    uh_idx   = ultraharmonic_indices(k);
    h_hi_idx = harmonic_indices(uh_upper_harm(k));

    bb_start = max(1,           uh_idx   + noise_exclusion + 1);
    bb_end   = min(length(pxx), h_hi_idx - noise_exclusion - 1);

    if bb_end <= bb_start
        warning('Empty BB region for U%.1f', uh_orders(k));
        continue
    end

    noise_psd = pxx(bb_start:bb_end);

    broadband_noise_power(k)    = mean(noise_psd * df);
    broadband_noise_power_dB(k) = 10*log10(broadband_noise_power(k));
    bb_bounds(k, :)              = [bb_start, bb_end];

    fprintf('U%.1f → H%d: [%.3f – %.3f] f0, BB = %.2f dB\n', ...
        uh_orders(k), uh_upper_harm(k)-1, ...
        f(bb_start)/drive_frequency, f(bb_end)/drive_frequency, ...
        broadband_noise_power_dB(k));
end

noise_floor_dB = mean(broadband_noise_power_dB);

%% === PLOT ===
f_norm      = f / drive_frequency;
spectrum_dB = 10*log10(pxx + eps);
dark_green  = [0 0.5 0];

figure('Color', 'w', 'Position', [200 200 900 600]);
ax = axes;
hold(ax, 'on');

% --- Full spectrum (gray) ---
plot(ax, f_norm, spectrum_dB, ...
    'Color', [0.75 0.75 0.75], 'LineWidth', 1.2, ...
    'DisplayName', 'Full Acoustic Spectrum');

% --- BB highlight: gap from each UH right edge to next harmonic ---
first_plotted = false;
for k = 1:length(uh_orders)
    bb_start = bb_bounds(k, 1);
    bb_end   = bb_bounds(k, 2);
    if bb_end <= bb_start, continue; end

    seg           = NaN(size(spectrum_dB));
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

% --- Noise floor line ---
yline(ax, noise_floor_dB, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Broadband Noise Floor = %.2f dB', noise_floor_dB));

% --- Harmonics (red circles) ---
plot(ax, f_norm(harmonic_indices), ...
    20*log10(harmonic_amplitudes + eps), ...
    'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, ...
    'DisplayName', 'Fundamental (f_0) and Harmonics');

% --- Ultraharmonics (green triangles) ---
valid_idx = ultraharmonic_indices > 0;
plot(ax, f_norm(ultraharmonic_indices(valid_idx)), ...
    20*log10(ultraharmonic_amplitudes(valid_idx) + eps), ...
    '^', 'Color', dark_green, 'MarkerFaceColor', dark_green, 'MarkerSize', 8, ...
    'DisplayName', 'Ultraharmonics');

% --- Labels ---
label_offset = 5;
for k = 1:n_harmonics
    if k == 1, label = 'f_0'; else, label = sprintf('H%d', k-1); end
    text(f_norm(harmonic_indices(k)), ...
        20*log10(harmonic_amplitudes(k)+eps) + label_offset, label, ...
        'Color', 'r', 'FontWeight', 'bold', 'FontSize', 12, ...
        'HorizontalAlignment', 'center');
end

for k = 1:length(uh_orders)
    if ultraharmonic_indices(k) <= 0, continue; end
    text(f_norm(ultraharmonic_indices(k)), ...
        20*log10(ultraharmonic_amplitudes(k)+eps) + label_offset, ...
        sprintf('U%.1f', uh_orders(k)), ...
        'Color', dark_green, 'FontWeight', 'bold', 'FontSize', 12, ...
        'HorizontalAlignment', 'center');
end

% --- Formatting ---
title(ax, 'Normalized Cavitation Spectrum Showing Harmonics, Ultraharmonics, and Broadband Noise', ...
    'FontSize', 14, 'FontWeight', 'bold');
xlabel(ax, 'Normalized Frequency  (f / f_0)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel(ax, 'Acoustic Power Spectral Density (dB)', 'FontSize', 13, 'FontWeight', 'bold');
set(ax, 'FontSize', 12, 'LineWidth', 1.2);
xlim(ax, [0 6]);
ylim(ax, [-150 -50]);
grid(ax, 'on');
box(ax, 'on');
legend(ax, 'Location', 'northeast');


