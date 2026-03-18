% in-vitro PCD Dataset (simplified)
clear; clc; close all;

% CHANGE THIS to your data folder
data_path = '/Users/shwethasuresh/Desktop/capstone/data/1_and_2_2026-01-16_single_tube_sweep/pcd/100x_01';   

fprintf('Using data folder: %s\n', data_path);

if ~exist(data_path, 'dir')
    error('Data folder not found: %s', data_path);
end

% Get all .mat files from the specified folder
mat_files = dir(fullfile(data_path, '*.mat'));

if isempty(mat_files)
    error('No .mat files found in: %s', data_path);
end

% Take every other file (alternating)
mat_files_alt = mat_files(2:2:end);

fprintf('\nFound %d alternating .mat files:\n', length(mat_files_alt));
for i = 1:length(mat_files_alt)
    fprintf('%d. %s\n', i, mat_files_alt(i).name);
end

file_number = 32;  % <-- change this as needed

if file_number > length(mat_files_alt)
    error('file_number exceeds available files.');
end

file_to_load = fullfile(data_path, mat_files_alt(file_number).name);
reader = load(file_to_load);

fprintf('\nData loaded successfully\n');

fprintf('\n=== DATA INFORMATION ===\n');
fprintf('File name: %s\n', mat_files_alt(file_number).name);
fprintf('Tstart: %d s\n', reader.Tstart);
fprintf('Tinterval: %d s\n', reader.Tinterval);
fprintf('Extra Samples: %.0f\n', reader.ExtraSamples);
fprintf('Requested Length: %.2f mm\n', reader.RequestedLength);
fprintf('Length: %.2f mm\n', reader.Length);
fprintf('Number of versions: %.0f\n', length(reader.Version));
fprintf('A length: %.0f\n', length(reader.A));
fprintf('B length: %.0f\n', length(reader.B));

%% look at raw B signal
%time t
t = reader.Tstart + (0:length(reader.B)-1) * reader.Tinterval;

%% apply fft to data
%pwelch formula
[pxx, f] = pwelch(reader.B, [], [], [], 1/reader.Tinterval);

% signal analysis
fprintf('\n=== SIGNAL ANALYSIS ===\n');
    
% signal statistics
fprintf('signal statistics:\n');
fprintf('  mean: %.4f mV\n', mean(reader.B));
fprintf('  standard deviation: %.4f mV\n', std(reader.B));
fprintf('  peak-to-peak: %.4f mV\n', max(reader.B) - min(reader.B));
fprintf('  RMS: %.4f mV\n', rms(reader.B));

%% find harmonics and ultraharmonics
% === Harmonic and Ultraharmonic Analysis using Welch Spectrum ===
WelchSpec = sqrt(pxx);

% Fundamental and harmonic setup
drive_frequency = 0.5e6;   % f0
n_harmonics = 4;           % f0, H1, H2, H3
harmonic_freqs = drive_frequency * (1:n_harmonics);

% Harmonic detection (Welch-based)
harmonic_indices = zeros(1, n_harmonics);
for k = 1:n_harmonics
    [~, idx] = min(abs(f - harmonic_freqs(k)));
    harmonic_indices(k) = idx;
end
harmonic_amplitudes = WelchSpec(harmonic_indices);

% Ultraharmonics: U1.5, U2.5, U3.5
uh_orders = [1.5 2.5 3.5];
ultraharmonic_freqs = uh_orders * drive_frequency;

ultraharmonic_indices = [];
ultraharmonic_amplitudes = [];

for k = 1:length(uh_orders)
    % Skip if frequency is above Welch range
    if ultraharmonic_freqs(k) > max(f)
        warning('Ultraharmonic U%.1f above Nyquist – skipped', uh_orders(k));
        continue
    end

    [~, idx] = min(abs(f - ultraharmonic_freqs(k)));

    ultraharmonic_indices(end+1) = idx; 
    ultraharmonic_amplitudes(end+1) = WelchSpec(idx);
end


%% Ultraharmonics bins
uh_orders = [1.5 2.5 3.5];      % ultraharmonic orders
drive_frequency = 0.5e6;        % fundamental frequency
ultraharmonic_freqs = uh_orders * drive_frequency;

% Set your search range in Hz (e.g., ±smth Hz around target)
df = f(2) - f(1);
% 1000 was arbitrary
bin_range = 1000; 
search_range = bin_range * df;  

% Preallocate arrays
ultraharmonic_indices = zeros(1, length(uh_orders));
ultraharmonic_amplitudes = zeros(1, length(uh_orders));

for k = 1:length(uh_orders)
    f_target = ultraharmonic_freqs(k);
    
    % Define search window
    idx_window = find(f >= f_target - search_range & f <= f_target + search_range);
    
    if isempty(idx_window)
        warning('No points found in search range for U%.1f', uh_orders(k));
        continue
    end
    
    % Find max within window
    [max_val, max_idx] = max(WelchSpec(idx_window));
    
    ultraharmonic_indices(k) = idx_window(max_idx);         % index of max
    ultraharmonic_amplitudes(k) = max_val;                 % amplitude at max
end

%% Print the results
fprintf('\n=== Ultraharmonics Peak Detection ===\n');
for k = 1:length(uh_orders)
    fprintf('U%.1f: Frequency = %.6f MHz, Amplitude = %.12f, Amplitude(dB) = %.2f dB\n', ...
        uh_orders(k), f(ultraharmonic_indices(k))/1e6, ...
        ultraharmonic_amplitudes(k), 20*log10(ultraharmonic_amplitudes(k)));
end
%% take the integral of harmonics & ultraharmonics to extract power which we want as the features
% --- Frequency resolution ---
df = f(2) - f(1);  % assuming uniform spacing

% --- Define a small bin range around each harmonic ---

fprintf('\n=== Total Harmonic Power (dB) ===\n');
for k = 1:n_harmonics
    if k == 1
        label = 'f0';
    else
        label = sprintf('H%d', k-1);
    end
    
    % Determine bin indices around the harmonic peak
    startBin = max(1, harmonic_indices(k) - bin_range);
    endBin   = min(length(pxx), harmonic_indices(k) + bin_range);
    
    % Integrate/sum the power under the harmonic
    harmonic_power = sum(pxx(startBin:endBin)) * df;
    
    % Convert to dB
    harmonic_power_dB = 10*log10(harmonic_power);  % use 10*log10 for power
    
    fprintf('%s: Frequency = %.3f MHz, Total power = %.14f mV^2, Total Power = %.2f dB\n', ...
            label, f(harmonic_indices(k))/1e6, harmonic_power, harmonic_power_dB);
end


fprintf('\n=== Ultraharmonic Total Power (dB) ===\n');
for k = 1:length(ultraharmonic_indices)
    % Determine bin indices around the ultraharmonic peak
    startBin = max(1, ultraharmonic_indices(k) - bin_range);
    endBin   = min(length(pxx), ultraharmonic_indices(k) + bin_range);
    
    % Integrate/sum the power under the ultraharmonic
    ultraharmonic_power = sum(pxx(startBin:endBin)) * df;
    
    % Convert to dB
    ultraharmonic_power_dB = 10*log10(ultraharmonic_power);
    
    fprintf('U%.1f: Frequency = %.3f MHz, Total Power = %.14f mv^2, Total Power = %.2f dB\n', ...
            uh_orders(k), f(ultraharmonic_indices(k))/1e6, ultraharmonic_power, ultraharmonic_power_dB);
end
%% Calculate Broadband Noise (restricted to 1.5f0 – 3f0)
drive_frequency = 0.5e6;       % fundamental frequency (Hz)
f_norm = f / drive_frequency;  % normalized frequency

% Define noise frequency range
noise_range_min = 1.5;  % 1.5 f0
noise_range_max = 3.0;  % 3 f0

df = f(2) - f(1);  % frequency resolution

% storage
broadband_noise_power = zeros(1,length(ultraharmonic_indices));
broadband_noise_power_dB =zeros(1,length(ultraharmonic_indices));

noise_exclusion = 50;  % bins excluded around UH peak

fprintf('\n=== Broadband Noise per Ultraharmonic ===\n');

for k = 1:length(ultraharmonic_indices)

    center_idx = ultraharmonic_indices(k);

    % Define window around ultraharmonic
    startBin = max(1, center_idx - bin_range);
    endBin   = min(length(pxx), center_idx + bin_range);

    % Initialize mask for this window
    local_mask = true(endBin-startBin+1,1);

    % Exclude the ultraharmonic peak region
    peak_start = max(startBin, center_idx - noise_exclusion);
    peak_end   = min(endBin, center_idx + noise_exclusion);

    local_mask((peak_start-startBin+1):(peak_end-startBin+1)) = false;

    % Extract PSD values
    local_psd = pxx(startBin:endBin);
    noise_psd = local_psd(local_mask);

    % Convert PSD → power
    noise_power_bins = noise_psd * df;

    % Average broadband noise
    broadband_noise_power(k) = mean(noise_power_bins);
    broadband_noise_power_dB(k) = 10*log10(broadband_noise_power(k));

    fprintf('U%.1f Broadband Noise: %.12f mV^2 (%.2f dB)\n', ...
        uh_orders(k), broadband_noise_power(k), broadband_noise_power_dB(k));

end
%% Normalized Harmonics + Ultraharmonics + Broadband Noise

drive_frequency = 0.5e6;        
f_norm = f / drive_frequency;   

% Avoid log of zero
spectrum_dB = 10*log10(pxx + eps);

% === Compute broadband noise floor (FIXED) ===
broadband_noise_avg_power_dB = mean(broadband_noise_power_dB);

% === Create figure ===
figure('Color','w','Position',[200 200 900 600]);
ax = axes;
hold(ax,'on');

%% === Plot full spectrum ===
plot(ax, f_norm, spectrum_dB, ...
     'Color', [0.75 0.75 0.75], ...
     'LineWidth', 1.2, ...
     'DisplayName','Full Acoustic Spectrum');

%% === Broadband noise region (1.5–3 f0) ===
noise_range_min = 1.5;  
noise_range_max = 3.0;  

freq_mask = (f_norm >= noise_range_min) & (f_norm <= noise_range_max);

% Reduce exclusion window (more reasonable)
noise_exclusion_range = 50;

peak_mask = true(size(pxx));
all_peak_indices = [harmonic_indices, ultraharmonic_indices];

for k = 1:length(all_peak_indices)
    idx = all_peak_indices(k);
    if idx <= 0, continue; end
    
    startBin = max(1, idx - noise_exclusion_range);
    endBin   = min(length(pxx), idx + noise_exclusion_range);
    peak_mask(startBin:endBin) = false;
end

final_noise_mask = freq_mask & peak_mask;

noise_spectrum = spectrum_dB;
noise_spectrum(~final_noise_mask) = NaN;

plot(ax, f_norm, noise_spectrum, ...
     'b', 'LineWidth',2, ...
     'DisplayName','Broadband Noise (1.5f_0 – 3f_0)');

%% === Noise floor line (FIXED) ===
yline(ax, broadband_noise_avg_power_dB, 'k--', ...
      'LineWidth',2, ...
      'DisplayName',sprintf('Noise Floor = %.2f dB', broadband_noise_avg_power_dB));

%% === Harmonics ===
plot(ax, f_norm(harmonic_indices), ...
     20*log10(harmonic_amplitudes + eps), ...
     'ro', ...
     'MarkerFaceColor','r', ...
     'MarkerSize',8, ...
     'DisplayName','Harmonics');

%% === Ultraharmonics (safe indexing) ===
valid_idx = ultraharmonic_indices > 0;

dark_green = [0 0.5 0];

plot(ax, f_norm(ultraharmonic_indices(valid_idx)), ...
     20*log10(ultraharmonic_amplitudes(valid_idx) + eps), ...
     '^', ...
     'Color', dark_green, ...
     'MarkerFaceColor', dark_green, ...
     'MarkerSize',8, ...
     'DisplayName','Ultraharmonics');

%% === Labels ===
label_offset = 5;
annotation_fontsize = 12;

% Harmonics labels
for k = 1:n_harmonics
    if k == 1
        label = 'f_0';
    else
        label = sprintf('H%d', k-1);
    end
    
    text(f_norm(harmonic_indices(k)), ...
         20*log10(harmonic_amplitudes(k)+eps) + label_offset, ...
         label, ...
         'Color','r', ...
         'FontWeight','bold', ...
         'FontSize',annotation_fontsize, ...
         'HorizontalAlignment','center');
end

% Ultraharmonic labels
for k = 1:length(uh_orders)
    if ultraharmonic_indices(k) <= 0, continue; end
    
    text(f_norm(ultraharmonic_indices(k)), ...
         20*log10(ultraharmonic_amplitudes(k)+eps) + label_offset, ...
         sprintf('U%.1f', uh_orders(k)), ...
         'Color', dark_green, ...
         'FontWeight','bold', ...
         'FontSize',annotation_fontsize, ...
         'HorizontalAlignment','center');
end

%% Normalized Harmonics + Ultraharmonics + Broadband Noise

drive_frequency = 0.5e6;        
f_norm = f / drive_frequency;   

% Avoid log of zero
spectrum_dB = 10*log10(pxx + eps);

% === Create figure ===
figure('Color','w','Position',[200 200 900 600]);
ax = axes;
hold(ax,'on');

%% === Plot full spectrum ===
plot(ax, f_norm, spectrum_dB, ...
     'Color', [0.75 0.75 0.75], ...
     'LineWidth', 1.2, ...
     'DisplayName','Full Acoustic Spectrum');

%% === Broadband noise region (1.5–3 f0) ===
noise_range_min = 1.5;  
noise_range_max = 3.0;  

freq_mask = (f_norm >= noise_range_min) & (f_norm <= noise_range_max);

% Reasonable exclusion window
noise_exclusion_range = 50;

peak_mask = true(size(pxx));
all_peak_indices = [harmonic_indices, ultraharmonic_indices];

for k = 1:length(all_peak_indices)
    idx = all_peak_indices(k);
    if idx <= 0, continue; end
    
    startBin = max(1, idx - noise_exclusion_range);
    endBin   = min(length(pxx), idx + noise_exclusion_range);
    peak_mask(startBin:endBin) = false;
end

final_noise_mask = freq_mask & peak_mask;

noise_spectrum = spectrum_dB;
noise_spectrum(~final_noise_mask) = NaN;

plot(ax, f_norm, noise_spectrum, ...
     'b', 'LineWidth',2, ...
     'DisplayName','Broadband Noise (1.5f_0 – 3f_0)');

%% === Harmonics ===
plot(ax, f_norm(harmonic_indices), ...
     20*log10(harmonic_amplitudes + eps), ...
     'ro', ...
     'MarkerFaceColor','r', ...
     'MarkerSize',8, ...
     'DisplayName','Harmonics');

%% === Ultraharmonics (safe indexing) ===
valid_idx = ultraharmonic_indices > 0;

dark_green = [0 0.5 0];

plot(ax, f_norm(ultraharmonic_indices(valid_idx)), ...
     20*log10(ultraharmonic_amplitudes(valid_idx) + eps), ...
     '^', ...
     'Color', dark_green, ...
     'MarkerFaceColor', dark_green, ...
     'MarkerSize',8, ...
     'DisplayName','Ultraharmonics');

%% === Labels ===
label_offset = 5;
annotation_fontsize = 12;

% Harmonics labels
for k = 1:n_harmonics
    if k == 1
        label = 'f_0';
    else
        label = sprintf('H%d', k-1);
    end
    
    text(f_norm(harmonic_indices(k)), ...
         20*log10(harmonic_amplitudes(k)+eps) + label_offset, ...
         label, ...
         'Color','r', ...
         'FontWeight','bold', ...
         'FontSize',annotation_fontsize, ...
         'HorizontalAlignment','center');
end

% Ultraharmonic labels
for k = 1:length(uh_orders)
    if ultraharmonic_indices(k) <= 0, continue; end
    
    text(f_norm(ultraharmonic_indices(k)), ...
         20*log10(ultraharmonic_amplitudes(k)+eps) + label_offset, ...
         sprintf('U%.1f', uh_orders(k)), ...
         'Color', dark_green, ...
         'FontWeight','bold', ...
         'FontSize',annotation_fontsize, ...
         'HorizontalAlignment','center');
end

%% === Formatting ===
title('Normalized Cavitation Spectrum', ...
      'FontSize',16,'FontWeight','bold')

xlabel('Normalized Frequency (f / f_0)', ...
       'FontSize',14,'FontWeight','bold')

ylabel('Power Spectral Density (dB)', ...
       'FontSize',14,'FontWeight','bold')

set(ax,'FontSize',12,'LineWidth',1.2)
xlim([0 6])
grid on
box on

legend('Location','northeast')
