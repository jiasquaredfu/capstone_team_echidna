%% Dataset 2 Analysis 
% Capstone II - Team E(Chidna)
% JiaJia Fu 
% 4/10/26

clear;
clc;
close all;

% Add path of imaging function 
function_path = '/Users/jiasquared/Desktop/CODING/capstone_team_echidna/04_Projects/';
addpath(function_path);

% Add path of Philips-ATL L12-5 Probe Data to read 
probe_path = ['/Users/jiasquared/Desktop/CODING/capstone_team_echidna' ...
    '/1_and_2_2026-01-16_single_tube_sweep/probe'];

clean_data = extract_probe_data(probe_path, 1);
fs = clean_data.fs;

%% Conversion to frequency domain for middle channel frame 1 example 
frame_1 = clean_data.rf_data(:, 64, 1);
[pxx, f] = pwelch(frame_1, [], [], [], fs);
pxx_db   = 10*log10(pxx);

%% Plotting 100x using Nicole's Thresholding Code from ambient PCD files
clear;
clc;
probe_path    = '/Users/jiasquared/Desktop/CODING/capstone_team_echidna/1_and_2_2026-01-16_single_tube_sweep/probe'; % want /PCI_100x_01.pacq
threshold_N = 1;  

fprintf('Using data folder:    %s\n', probe_path);

probe_files = dir(fullfile(probe_path, '*.pacq*'));

fprintf('Found probe files:\n');
for i = 1:length(probe_files)
    fprintf('  %d. %s\n', i, probe_files(i).name);
end

target_name = '100x_01';
ambient_name = 'ambient';
probe_idx   = find(contains({probe_files.name}, target_name), 1);
ambient_idx   = find(contains({probe_files.name}, ambient_name), 1);

if isempty(probe_idx)
    error('Could not find probe file matching "%s"', target_name);
end
fprintf('\nUsing probe file: %s\n', probe_files(probe_idx).name);
fprintf('\nUsing probe ambient file: %s\n', probe_files(ambient_idx).name);

clean_data = extract_probe_data(probe_path, probe_idx); %extract 100x_1 file for comparison 
fs_32 = clean_data.fs;

fprintf('Probe fs:      %.0f Hz\n',  fs_32);
fprintf('Probe frames:  %d\n',       size(clean_data.rf_data, 3));
fprintf('Probe channels: %d\n',      size(clean_data.rf_data, 2));

ambient_data = extract_probe_data(probe_path, ambient_idx); %extract ambient file for comparison 
a_fs = ambient_data.fs;


fprintf('Ambient fs:      %.0f Hz\n',  a_fs);
fprintf('Ambient frames:  %d\n',       size(ambient_data.rf_data, 3));
fprintf('Ambient channels: %d\n',      size(ambient_data.rf_data, 2));



%% Load probe and ambient data — middle channel, frame 32

frame_32  = clean_data.rf_data(:, 64, 1); % highest pressure at 1, best plot 
[pxx_32, f_32] = pwelch(frame_32, [], [], [], fs_32);
WelchSpec = pxx_32; %just in case don't break code downstream
pxx_db_32 = 10*log10(pxx_32);
df_32 = f_32(2) - f_32(1);
f0 = 0.5e6;
f_norm_32 = f_32 / f0;
bin_range = 1000;

% === SIGNAL STATS ===
fprintf('\n=== SIGNAL ANALYSIS ===\n');
fprintf('  Mean:          %.4f mV\n', mean(frame_32));
fprintf('  Std deviation: %.4f mV\n', std(frame_32));
fprintf('  Peak-to-peak:  %.4f mV\n', max(frame_32) - min(frame_32))

ambient_frame = ambient_data.rf_data(:, 64, 1);  
[pxx_amb, f_amb] = pwelch(ambient_frame, [], [], [], a_fs);
amb_spec = pxx_amb;  %just in case don't break code downstream
amb_df   = f_amb(2) - f_amb(1);
pxx_amb_db = 10*log10(pxx_amb);
amb_sr   = bin_range * amb_df;

fprintf('\n=== AMBIENT SIGNAL STATS ===\n');
fprintf('  Mean:          %.4f mV\n', mean(ambient_frame));
fprintf('  Std deviation: %.4f mV\n', std(ambient_frame));
fprintf('  Peak-to-peak:  %.4f mV\n', max(ambient_frame) - min(ambient_frame));

%% Harmonic Ultraharmonic Detection Test

% === HARMONIC DETECTION PARAMS ===
harmonics =  6:9; 
n_harmonics = length(harmonics); 
harmonic_freqs = f0 * harmonics;

uh_orders           = [6.5 7.5 8.5];
ultraharmonic_freqs = uh_orders * f0;
search_range        = bin_range * df_32;

harmonic_indices = zeros(1, n_harmonics);
for k = 1:n_harmonics
    [~, idx] = min(abs(f_32 - harmonic_freqs(k)));
    harmonic_indices(k) = idx;
end
harmonic_amplitudes = WelchSpec(harmonic_indices);

ultraharmonic_indices    = zeros(1, length(uh_orders));
ultraharmonic_amplitudes = zeros(1, length(uh_orders));

for k = 1:length(uh_orders)
    f_target   = ultraharmonic_freqs(k);
    idx_window = find(f_32 >= f_target - search_range & f_32 <= f_target + search_range);
    if isempty(idx_window)
        warning('No points found in search range for U%.1f', uh_orders(k));
        continue
    end
    [max_val, max_idx]          = max(WelchSpec(idx_window));
    ultraharmonic_indices(k)    = idx_window(max_idx);
    ultraharmonic_amplitudes(k) = max_val;
end

disp(harmonic_amplitudes)
disp(ultraharmonic_indices)
disp(ultraharmonic_amplitudes)

%% H and UH Thresholding
% === AMBIENT H AND UH DETECTION ===

% UH peak amplitudes (dB) from the matched ambient file
amb_uh_peak_dB = NaN(1, length(uh_orders));
for k = 1:length(uh_orders)
    f_target   = uh_orders(k) * f0;
    idx_window = find(f_amb >= f_target - amb_sr & f_amb <= f_target + amb_sr);
    if isempty(idx_window)
        warning('No ambient points found for U%.1f', uh_orders(k));
        continue
    end
    peak_amp           = max(amb_spec(idx_window));
    amb_uh_peak_dB(k)  = 10*log10(peak_amp);
end

amb_uh_std_dB = NaN(1, length(uh_orders));
for k = 1:length(uh_orders)
    f_target   = uh_orders(k) * f0;
    idx_window = find(f_amb >= f_target - amb_sr & f_amb <= f_target + amb_sr);
    if isempty(idx_window), continue; end
    window_vals_dB      = 10*log10(amb_spec(idx_window) + eps);
    amb_uh_std_dB(k)    = std(window_vals_dB);
end

threshold_dB = amb_uh_peak_dB + threshold_N * amb_uh_std_dB;

% Target file UH peak dB (for comparison)
target_uh_dB = 10*log10(ultraharmonic_amplitudes);

% Exceeds threshold flag (UH)
exceeds = target_uh_dB > threshold_dB;


%% === HARMONIC THRESHOLD COMPUTATION ===

% Use gaps between harmonics, excluding ±bin_range bins around each harmonic
amb_noise_mask = true(size(f_amb));
for k = 1:n_harmonics
    [~, h_idx] = min(abs(f_amb - harmonic_freqs(k)));
    lo = max(1, h_idx - bin_range);
    hi = min(length(f_amb), h_idx + bin_range);
    amb_noise_mask(lo:hi) = false;
end
% Also exclude UH regions
for k = 1:length(uh_orders)
    f_target = uh_orders(k) * f0;
    [~, u_idx] = min(abs(f_amb - f_target));
    lo = max(1, u_idx - bin_range);
    hi = min(length(f_amb), u_idx + bin_range);
    amb_noise_mask(lo:hi) = false;
end

amb_noise_vals_dB  = 10*log10(amb_spec(amb_noise_mask));
amb_noise_mean_dB  = mean(amb_noise_vals_dB);
amb_noise_std_dB   = std(amb_noise_vals_dB);
harmonic_threshold_dB_scalar = amb_noise_mean_dB + threshold_N * amb_noise_std_dB;

harmonic_threshold_dB = repmat(harmonic_threshold_dB_scalar, 1, n_harmonics);
amb_harm_peak_dB      = repmat(amb_noise_mean_dB, 1, n_harmonics);
amb_harm_std_dB       = repmat(amb_noise_std_dB,  1, n_harmonics);

target_harm_dB = 10*log10(harmonic_amplitudes);
exceeds_harm = false(1, n_harmonics);
for k = 1:n_harmonics
    exceeds_harm(k) = target_harm_dB(k) > harmonic_threshold_dB(k);
end

%% === PRINT THRESHOLD TABLES ===

% --- Harmonic table (H2, H3, H4) ---
fprintf('\n=== Harmonic Threshold Summary ===\n');
fprintf('Threshold = Ambient Noise Floor Mean + %d × Std  (inter-harmonic gaps)\n', threshold_N);
fprintf('Ambient noise floor: Mean = %.2f dB,  Std = %.2f dB,  Threshold = %.2f dB\n\n', ...
    amb_noise_mean_dB, amb_noise_std_dB, harmonic_threshold_dB_scalar);
fprintf('%-6s  %12s  %10s  %8s\n', 'H', 'Target Peak', 'Threshold', 'Exceeds?');
fprintf('%s\n', repmat('-', 1, 45));
for k = 1:n_harmonics
    flag = 'NO';
    if exceeds_harm(k), flag = 'YES ***'; end
    fprintf('H%-5d  %+12.2f  %+10.2f  %s\n', ...
        harmonics(k), target_harm_dB(k), harmonic_threshold_dB(k), flag);
end
fprintf('%s\n', repmat('-', 1, 45));

% --- Ultraharmonic table ---
fprintf('\n=== Ultraharmonic Threshold Summary ===\n');
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
fprintf('All values in dB (10·log10 of Welch amplitude spectrum)\n');


%% === SAFETY ASSESSMENT: U6.5 Broadband Region ===
% Bandpass filter both target and ambient to the U6.5 broadband region
% [6.5*f0 + bin_range*df  →  7.0*f0 - bin_range*df]

f0 = 0.5e6;          % fundamental frequency
bin_range = 1000;    % bin range for summing


% --- Compute bandpass limits ---
bb_lo_raw = 6.5*f0 + bin_range*df_32;
bb_hi_raw = 7.0*f0 - bin_range*df_32;

% Ensure bb_lo < bb_hi
bb_lo_Hz = min(bb_lo_raw, bb_hi_raw);
bb_hi_Hz = max(bb_lo_raw, bb_hi_raw);

% Normalize to Nyquist
bb_lo = bb_lo_Hz / (fs_32/2);
bb_hi = bb_hi_Hz / (fs_32/2);

% Clamp to valid range
bb_lo = max(bb_lo, 1e-4);
bb_hi = min(bb_hi, 1-1e-4);

% --- Design 4th-order Butterworth bandpass filter ---
[b_bp, a_bp] = butter(4, [bb_lo, bb_hi], 'bandpass');

% --- Apply filter ---
target_filtered  = filtfilt(b_bp, a_bp, frame_32);
ambient_filtered = filtfilt(b_bp, a_bp, ambient_frame);

% --- Compute standard deviations ---
target_bb_std  = std(target_filtered);
ambient_bb_std = std(ambient_filtered);

% --- Thresholds ---
thresh_1s = ambient_bb_std * 2;
thresh_2s = ambient_bb_std * 3;
thresh_3s = ambient_bb_std * 4;

% --- Determine safety level ---
if target_bb_std > thresh_3s
    safety_label = 'UNSAFE (> 3σ)';
elseif target_bb_std > thresh_2s
    safety_label = 'UNSAFE (> 2σ)';
elseif target_bb_std > thresh_1s
    safety_label = 'UNSAFE (> 1σ)';
else
    safety_label = 'SAFE';
end

% --- Display results ---
fprintf('\n=== Safety Assessment: U6.5 Broadband Region ===\n');
fprintf('Bandpass range:        %.4f – %.4f MHz\n', bb_lo_Hz/1e6, bb_hi_Hz/1e6);
fprintf('Ambient BB std:        %.6f\n', ambient_bb_std);
fprintf('Target  BB std:        %.6f\n', target_bb_std);

fprintf('\n  Threshold 1σ:  %.6f  →  %s\n', thresh_1s, string_exceed(target_bb_std > thresh_1s));
fprintf('  Threshold 2σ:  %.6f  →  %s\n', thresh_2s, string_exceed(target_bb_std > thresh_2s));
fprintf('  Threshold 3σ:  %.6f  →  %s\n', thresh_3s, string_exceed(target_bb_std > thresh_3s));

fprintf('\n  *** OVERALL: %s ***\n', safety_label);
fprintf('%s\n', repmat('-', 1, 50));

% === Nested helper function ===
function s = string_exceed(cond)
    if cond
        s = 'exceeded';
    else
        s = 'within';
    end
end

%% Plot with all labels

figure; grid on; hold on;
dark_green  = [0 0.5 0];

% --- Full spectrum ---
plot(f_norm_32, pxx_db_32, 'Color',[0.75 0.75 0.75], 'LineWidth',1.2, ...
     'DisplayName','Full Acoustic Spectrum');

% --- Broadband noise regions ---
shrink = 0.05;  % small margin instead of bin_range_norm
bb_regions = [6.5 + shrink, 7.0 - shrink;
              7.5 + shrink, 8.0 - shrink;
              8.5 + shrink, 9.0 - shrink];

first_plotted = false;
for k = 1:size(bb_regions,1)
    mask = f_norm_32 >= bb_regions(k,1) & f_norm_32 <= bb_regions(k,2);
    if ~any(mask), continue; end
    seg = NaN(size(pxx_db_32));
    seg(mask) = pxx_db_32(mask);
    if ~first_plotted
        plot(f_norm_32, seg, 'b', 'LineWidth', 2, 'DisplayName', 'Broadband Noise Regions');
        first_plotted = true;
    else
        plot(f_norm_32, seg, 'b', 'LineWidth', 2, 'HandleVisibility','off');
    end
end

% --- Harmonics (red circles) with labels ---
label_offset = 1.5;
for k = 1:n_harmonics
    idx = harmonic_indices(k);
    if idx <= 0, continue; end
    plot(f_norm_32(idx), 10*log10(harmonic_amplitudes(k)), ...
         'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 18, 'HandleVisibility','off');
    text(f_norm_32(idx), 10*log10(harmonic_amplitudes(k))+label_offset, ...
         sprintf('H%d', harmonics(k)), 'Color','r', 'FontWeight','bold', 'FontSize',18, ...
         'HorizontalAlignment','center');
end
% Add legend for harmonics
plot(NaN, NaN, 'ro', 'MarkerFaceColor','r', 'MarkerSize',18, 'DisplayName','Harmonics');

% --- Ultraharmonics (green triangles) with labels ---
for k = 1:length(uh_orders)
    idx = ultraharmonic_indices(k);
    if idx <= 0, continue; end
    plot(f_norm_32(idx), 10*log10(ultraharmonic_amplitudes(k)), ...
         '^', 'Color', dark_green, 'MarkerFaceColor', dark_green, 'MarkerSize', 18, 'HandleVisibility','off');
    text(f_norm_32(idx), 10*log10(ultraharmonic_amplitudes(k))+label_offset, ...
         sprintf('U%.1f', uh_orders(k)), 'Color', dark_green, 'FontWeight','bold', 'FontSize',18, ...
         'HorizontalAlignment','center');
end
% Add legend for ultraharmonics
plot(NaN, NaN, '^', 'Color', dark_green, 'MarkerFaceColor', dark_green, 'MarkerSize',18, 'DisplayName','Ultraharmonics');

% --- Harmonic threshold lines (H2, H3, H4) ---
harm_thr_colors = [0.85 0.3 0; 0.75 0.1 0; 0.60 0.0 0];  % colors
harm_labels = arrayfun(@(h) sprintf('H%d', h), harmonics, 'UniformOutput', false);
line_frac = 0.4;
for k = 2:n_harmonics
    if isnan(harmonic_threshold_dB(k)), continue; end
    f_ctr = harmonic_freqs(k)/f0;
    half_w = search_range/f0 * line_frac;
    x_seg = [f_ctr-half_w, f_ctr+half_w];
    y_seg = [harmonic_threshold_dB(k), harmonic_threshold_dB(k)];
    plot(x_seg, y_seg, '--', 'Color', harm_thr_colors(min(k-1,3),:), 'LineWidth',2.5, ...
         'DisplayName', sprintf('Threshold %s = %.1f dB', harm_labels{k}, harmonic_threshold_dB(k)));
    text(f_ctr, harmonic_threshold_dB(k)+2, sprintf('Thr %s', harm_labels{k}), ...
         'Color', harm_thr_colors(min(k-1,3),:), 'FontSize',18, 'HorizontalAlignment','center','FontWeight','bold');
end

% --- Ultraharmonic threshold lines ---
uh_colors = [0.85 0 0.85; 0.60 0 0.80; 0.30 0 0.70];
for k = 1:length(uh_orders)
    f_ctr = uh_orders(k);
    half_w = search_range/f0* line_frac;
    x_seg = [f_ctr-half_w, f_ctr+half_w];
    y_seg = [threshold_dB(k), threshold_dB(k)];
    plot(x_seg, y_seg, '--', 'Color', uh_colors(min(k,3),:), 'LineWidth',2.5, ...
         'DisplayName', sprintf('Threshold U%.1f = %.1f dB', uh_orders(k), threshold_dB(k)));
    text(f_ctr, threshold_dB(k)+2, sprintf('Thr U%.1f', uh_orders(k)), ...
         'Color', uh_colors(min(k,3),:), 'FontSize',18, 'HorizontalAlignment','center','FontWeight','bold');
end

% --- Formatting ---
title('Normalized Cavitation Spectrum — Harmonics, Ultraharmonics, Broadband Noise, and Thresholds', ...
      'FontSize',24,'FontWeight','bold');
xlabel('Normalized Frequency  (f / f_0)','FontSize',24,'FontWeight','bold');
ylabel('Acoustic Power Spectral Density (dB)','FontSize',24,'FontWeight','bold');
ylim([-65 -40]); xlim([5 10]);
grid on; box on;
legend('Location','northeast','FontSize',24);




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

    % --- Harmonics ---
    for i = 1:length(harmonics)
        target_f  = harmonics(i) * f0;
        [~, cidx] = min(abs(f - target_f));
        idx_range = max(1, cidx-search_window) : min(length(f), cidx+search_window);

        [~, local] = max(pxx_db(idx_range));
        peak_idx   = idx_range(local);

        h = harmonics(i);
        features.(sprintf('H%d_freq',  h)) = f(peak_idx);
        features.(sprintf('H%d_db',    h)) = pxx_db(peak_idx);
        features.(sprintf('H%d_power', h)) = trapz(f(idx_range), pxx(idx_range));
        harm_freq_vals(i) = f(peak_idx);
    end

    % --- Ultraharmonics ---
    for i = 1:length(ultraharmonics)
        target_f  = ultraharmonics(i) * f0;
        [~, cidx] = min(abs(f - target_f));
        idx_range = max(1, cidx-search_window) : min(length(f), cidx+search_window);

        [~, local] = max(pxx_db(idx_range));
        peak_idx   = idx_range(local);

        u_int = ultraharmonics(i) * 10;          % 6.5->65, stored as U65_freq etc.
        features.(sprintf('U%d_freq',  u_int)) = f(peak_idx);
        features.(sprintf('U%d_db',    u_int)) = pxx_db(peak_idx);
        features.(sprintf('U%d_power', u_int)) = trapz(f(idx_range), pxx(idx_range));
        ultra_freq_vals(i) = f(peak_idx);
    end

    % --- Broadband noise ---
    all_peak_freqs     = [harm_freq_vals(:); ultra_freq_vals(:)];
    noise_mask         = (f >= 9*f0) & (f <= 10.5*f0);
    for k = 1:length(all_peak_freqs)
        noise_mask(abs(f - all_peak_freqs(k)) <= 0.08e6) = false;
    end

    noise_psd                      = pxx(noise_mask);
    features.broadband_noise_db    = 10*log10(mean(noise_psd));
    features.broadband_noise_power = trapz(f(noise_mask), noise_psd);
end


function [thr1, thr2, thr3] = compute_pcd_threshold(ambient_folder, f0)
    bb_vals   = [];
    mat_files = dir(fullfile(ambient_folder, '*.mat'));

    if isempty(mat_files)
        warning('No ambient PCD .mat files found in %s', ambient_folder);
        thr1 = NaN; thr2 = NaN; thr3 = NaN;
        return
    end

    for i = 1:numel(mat_files)
        fullpath     = fullfile(mat_files(i).folder, mat_files(i).name);
        ambient_data = load(fullpath);

        % Known field layout: A = signal, Tinterval = sample period (s)
        if ~isfield(ambient_data, 'A') || ~isfield(ambient_data, 'Tinterval')
            warning('Unexpected fields in %s — skipping', mat_files(i).name);
            continue
        end

        rf_data = ambient_data.A;
        fs_val  = 1 / ambient_data.Tinterval;
        fprintf('File: %s | fs = %.0f Hz | size(A) = %s\n', ...
                mat_files(i).name, fs_val, mat2str(size(rf_data)));

        % A is a single column vector — one continuous signal, no frames
        [pxx, f] = pwelch(rf_data(:), [], [], [], fs_val);
        feat     = get_features(f, pxx, f0);
        bb_vals(end+1) = feat.broadband_noise_db; %#ok<AGROW>
    end

    bb_vals = bb_vals(~isnan(bb_vals));

    if isempty(bb_vals)
        warning('No broadband noise values extracted.');
        thr1 = NaN; thr2 = NaN; thr3 = NaN;
    else
        mu    = mean(bb_vals);
        sigma = std(bb_vals);
        thr1  = mu +   sigma;
        thr2  = mu + 2*sigma;
        thr3  = mu + 3*sigma;
        fprintf('Thresholds — k1: %.2f | k2: %.2f | k3: %.2f\n', thr1, thr2, thr3);
    end
end
