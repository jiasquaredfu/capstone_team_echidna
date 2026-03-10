% in-vitro PCD Dataset
clear; clc; close all;

fprintf('current folder: %s\n', pwd)

% add repo to pat
repo_path = '04_Projects/17_FUS_Instruments_Interface';
if exist(repo_path, 'dir')
    addpath(genpath(repo_path));
    fprintf('added repository to path\n');
else
    error('repo not found at: %s', repo_path);
end

% find .mat files in a different folder
mat_files = dir('*.mat');

% alternate files b/c they increase by pressure (every other one)
mat_files_alt = mat_files(1:2:end);

% display number and what matfiles are in mat_files_alt
fprintf('\nfound %d .mat files (alternating):\n', length(mat_files_alt));
for i = 1:length(mat_files_alt)
    fprintf('%d. %s\n', i, mat_files_alt(i).name);
end

% read file
file_number = 32;
reader = load(mat_files_alt(file_number).name);

% data is now stored in the reader object's properties
fprintf('\ndata loaded successfully into reader object\n');

% display basic information from the reader
fprintf('\n=== DATA INFORMATION ===\n');
fprintf('file name: %s\n', mat_files_alt(file_number).name);
fprintf('Tstart: %d s\n', reader.Tstart);
fprintf('Tinterval: %d s\n', reader.Tinterval);
fprintf('Number of Extra Samples: %.0f \n', reader.ExtraSamples);
fprintf('Requested Length: %.2f mm\n', reader.RequestedLength);
fprintf('Length: %.2f mm\n', reader.Length);
fprintf('Number of versions: %.2f \n', length(reader.Version));
fprintf('A length of table: %.2f \n', length(reader.A));
fprintf('B length of table: %.2f \n', length(reader.B));


%% look at raw B signal
%time t
t = reader.Tstart + (0:length(reader.B)-1) * reader.Tinterval;

% %raw signal just to see
% fig = figure('Name','Raw Signals','NumberTitle','off','Position',[100 100 1000 600]);
% tgroup = uitabgroup(fig);
% 
% % --- Tab 1: A signal ---
% tab1 = uitab(tgroup,'Title','A signal');
% ax1 = axes('Parent', tab1);
% plot(ax1, t, reader.A, 'b');
% title(ax1,'Raw A Signal')
% xlabel(ax1,'Time (s)')
% ylabel(ax1,'Amplitude')
% grid(ax1,'on')
% 
% % --- Tab 2: B signal ---
% tab2 = uitab(tgroup,'Title','B signal');
% ax2 = axes('Parent', tab2);
% plot(ax2, t, reader.B, 'r');
% title(ax2,'Raw B Signal')
% xlabel(ax2,'Time (s)')
% ylabel(ax2,'Amplitude')
% grid(ax2,'on')


%% apply fft to data
%pwelch formula
[pxx, f] = pwelch(reader.B, [], [], [], 1/reader.Tinterval);

%plot fft
fig = figure('Name','FFT B Signal','NumberTitle','off','Position',[100 100 1000 600]);
tgroup = uitabgroup(fig);

tab1 = uitab(tgroup,'Title','power spectrum');
ax1 = axes('Parent', tab1);
plot(f/1e6, 10*log10(pxx), 'B');
title(ax1,'pwelch applied to B signal')
xlabel(ax1,'frequencuy (MHz)')
ylabel(ax1,'power(dB)')
xlim([0 3])
grid(ax1,'on')


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

    ultraharmonic_indices(end+1) = idx; %#ok<SAGROW>
    ultraharmonic_amplitudes(end+1) = WelchSpec(idx);
end

% === Plot Welch Spectrum with Harmonics ===
tab2 = uitab(tgroup,'Title','Harmonics labeled');
ax2 = axes('Parent', tab2);

plot(ax2, f/1e6, 20*log10(WelchSpec));
hold(ax2,'on')

% Harmonic markers
plot(ax2, f(harmonic_indices)/1e6, ...
     20*log10(harmonic_amplitudes), ...
     'ro', 'MarkerFaceColor','r');

% Ultraharmonic markers
plot(ax2, f(ultraharmonic_indices)/1e6, ...
     20*log10(ultraharmonic_amplitudes), ...
     'go', 'MarkerFaceColor','g');

% Labels
label_offset = 0.05 * max(20*log10(WelchSpec));

for k = 1:n_harmonics
    if k == 1
        label = 'f0';
    else
        label = sprintf('H%d', k-1);
    end
    text(f(harmonic_indices(k))/1e6, ...
         20*log10(harmonic_amplitudes(k)) + label_offset, ...
         label, 'Color','r', 'FontWeight','bold', ...
         'HorizontalAlignment','center');
end

for k = 1:length(uh_orders)
    text(f(ultraharmonic_indices(k))/1e6, ...
         20*log10(ultraharmonic_amplitudes(k)) + label_offset, ...
         sprintf('U%.1f', uh_orders(k)), ...
         'Color','g', 'FontWeight','bold', ...
         'HorizontalAlignment','center');
end

% Axes formatting
title(ax2,'Welch Spectrum of B Signal with Harmonics')
xlabel(ax2,'Frequency (MHz)')
ylabel(ax2,'Magnitude (dB)')
xlim(ax2,[0 3])
grid(ax2,'on')

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
% === Plot PWELCH Spectrum with Harmonics and Ultraharmonics (markers only) ===
tab3 = uitab(tgroup,'Title','Ultraharmonics with U bins');
ax3 = axes('Parent', tab3);

% Plot full PWELCH spectrum (thin line)
plot(ax3, f/1e6, 20*log10(WelchSpec), 'b', 'LineWidth', 1, 'DisplayName','PWELCH');
hold(ax3,'on');

% Plot harmonics as markers only (red circles)
plot(ax3, f(harmonic_indices)/1e6, 20*log10(harmonic_amplitudes), 'ro', ...
     'MarkerFaceColor','r', 'MarkerSize',6, 'DisplayName','Harmonics');

% Plot ultraharmonics as markers only (green triangles)
plot(ax3, f(ultraharmonic_indices)/1e6, 20*log10(ultraharmonic_amplitudes), 'g^', ...
     'MarkerFaceColor','g', 'MarkerSize',6, 'DisplayName','Ultraharmonics');

% Add labels above peaks
label_offset = 0.05 * max(20*log10(WelchSpec));

for k = 1:n_harmonics
    if k == 1
        label = 'f0'; % fundamental
    else
        label = sprintf('H%d', k-1);
    end
    text(f(harmonic_indices(k))/1e6, ...
         20*log10(harmonic_amplitudes(k)) + label_offset, ...
         label, 'Color','r', 'FontWeight','bold', ...
         'HorizontalAlignment','center');
end

fprintf('\n=== Harmonic Values (dB) ===\n');
for k = 1:n_harmonics
    if k == 1
        label = 'f0';
    else
        label = sprintf('H%d', k-1);
    end
    % Convert amplitude to dB
    amplitude_dB = 20*log10(harmonic_amplitudes(k));
    
    fprintf('%s: Frequency = %.3f MHz, Amplitude = %.2f dB\n', ...
            label, f(harmonic_indices(k))/1e6, amplitude_dB);
end


for k = 1:length(ultraharmonic_indices)
    text(f(ultraharmonic_indices(k))/1e6, ...
         20*log10(ultraharmonic_amplitudes(k)) + label_offset, ...
         sprintf('U%.1f', uh_orders(k)), 'Color','g', 'FontWeight','bold', ...
         'HorizontalAlignment','center');
end

% Axes formatting
title(ax3,'PWELCH + Harmonics + Ultraharmonics')
xlabel(ax3,'Frequency (MHz)')
ylabel(ax3,'Magnitude (dB)')
xlim(ax3,[0 3])
grid(ax3,'on')
legend(ax3,'show')

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

%% normalization graph (f/f0) 
drive_frequency = 0.5e6;   % fundamental frequency (Hz)
f_norm = f / drive_frequency;

tab4 = uitab(tgroup,'Title','Normalized PWELCH');
ax4 = axes('Parent', tab4);

% Plot full PWELCH spectrum (thin line)
plot(ax4, f_norm, 20*log10(WelchSpec), 'b', ...
     'LineWidth', 1, 'DisplayName','PWELCH');
hold(ax4,'on');

% Plot harmonics as markers only (red circles)
plot(ax4, f_norm(harmonic_indices), ...
     20*log10(harmonic_amplitudes), 'ro', ...
     'MarkerFaceColor','r', ...
     'MarkerSize',6, ...
     'DisplayName','Harmonics');

% Plot ultraharmonics as markers only (green triangles)
plot(ax4, f_norm(ultraharmonic_indices), ...
     20*log10(ultraharmonic_amplitudes), 'g^', ...
     'MarkerFaceColor','g', ...
     'MarkerSize',6, ...
     'DisplayName','Ultraharmonics');

% Add labels above peaks
label_offset = 0.05 * max(20*log10(WelchSpec));

for k = 1:n_harmonics
    if k == 1
        label = 'f0';
    else
        label = sprintf('H%d', k-1);
    end
    
    text(f_norm(harmonic_indices(k)), ...
         20*log10(harmonic_amplitudes(k)) + label_offset, ...
         label, ...
         'Color','r', ...
         'FontWeight','bold', ...
         'HorizontalAlignment','center');
end

for k = 1:length(ultraharmonic_indices)
    text(f_norm(ultraharmonic_indices(k)), ...
         20*log10(ultraharmonic_amplitudes(k)) + label_offset, ...
         sprintf('U%.1f', uh_orders(k)), ...
         'Color','g', ...
         'FontWeight','bold', ...
         'HorizontalAlignment','center');
end

% Axes formatting
title(ax4,'PWELCH + Harmonics + Ultraharmonics (Normalized to f_0)')
xlabel(ax4,'Frequency / f_0 (f_0 = 0.5 MHz)')
ylabel(ax4,'Magnitude (dB)')
xlim(ax4,[0 6])   % adjust if needed
grid(ax4,'on')
legend(ax4,'show')

%% Calculate Broadband Noise (restricted to 1.5f0 – 3f0)
drive_frequency = 0.5e6;       % fundamental frequency (Hz)
f_norm = f / drive_frequency;  % normalized frequency

% Define noise frequency range
noise_range_min = 1.5;  % 1.5 f0
noise_range_max = 3.0;  % 3 f0

% Exclude regions around harmonics and ultraharmonics
noise_exclusion_range = 2000;  % bins to exclude around each peak

% Create mask for frequency bins within 1.5f0 to 3f0
freq_mask = (f_norm >= noise_range_min) & (f_norm <= noise_range_max);

% Create mask for all frequency bins initially
noise_mask = true(size(pxx));

% Combine all peak indices
all_peak_indices = [harmonic_indices, ultraharmonic_indices];

% Exclude regions around each peak
for k = 1:length(all_peak_indices)
    startBin = max(1, all_peak_indices(k) - noise_exclusion_range);
    endBin   = min(length(pxx), all_peak_indices(k) + noise_exclusion_range);
    noise_mask(startBin:endBin) = false;
end

% Combine with frequency mask
final_noise_mask = noise_mask & freq_mask;

% Get noise bins (power spectral density values)
noise_bins_psd = pxx(final_noise_mask);

% Calculate power of each noise bin
df = f(2) - f(1);  % frequency resolution
noise_bins_power = noise_bins_psd * df;  % convert PSD to power

% Take the average of the power values
broadband_noise_avg_power = mean(noise_bins_power);
broadband_noise_avg_power_dB = 10*log10(broadband_noise_avg_power);

fprintf('\n=== Broadband Noise (1.5f0 – 3f0) ===\n');
fprintf('Mean Power per Bin: %.12f mV^2\n', broadband_noise_avg_power);
fprintf('Mean Power per Bin (dB): %.2f dB\n', broadband_noise_avg_power_dB);
fprintf('Number of bins used for noise: %d / %d (%.1f%%)\n', ...
        sum(final_noise_mask), length(pxx), 100*sum(final_noise_mask)/length(pxx));

%% Normalized Harmonics + Ultraharmonics + Broadband Noise (Standalone Figure)
drive_frequency = 0.5e6;        
f_norm = f / drive_frequency;   
spectrum_dB = 10*log10(pxx);

% --------------------------------------------------
% Create New Figure
% --------------------------------------------------
figure('Color','w','Position',[200 200 900 600]);
ax = axes;
hold(ax,'on');

% --------------------------------------------------
% 1. Plot Full Spectrum
% --------------------------------------------------
plot(ax, f_norm, spectrum_dB, ...
     'Color', [0.75 0.75 0.75], ...
     'LineWidth', 1.2, ...
     'DisplayName','Full Acoustic Spectrum');

% --------------------------------------------------
% 2. Highlight Broadband Noise (1.5f0 – 3f0, excluding peaks)
% --------------------------------------------------
% Define normalized frequency range
noise_range_min = 1.5;  
noise_range_max = 3.0;  

% Mask for frequency range
freq_mask = (f_norm >= noise_range_min) & (f_norm <= noise_range_max);

% Mask for excluded peaks around harmonics and ultraharmonics
noise_exclusion_range = 2000;  % bins around peaks
peak_mask = true(size(pxx));
all_peak_indices = [harmonic_indices, ultraharmonic_indices];

for k = 1:length(all_peak_indices)
    startBin = max(1, all_peak_indices(k) - noise_exclusion_range);
    endBin   = min(length(pxx), all_peak_indices(k) + noise_exclusion_range);
    peak_mask(startBin:endBin) = false;
end

% Final noise mask
final_noise_mask = freq_mask & peak_mask;

% Prepare noise spectrum for plotting
noise_spectrum = spectrum_dB;
noise_spectrum(~final_noise_mask) = NaN;

plot(ax, f_norm, noise_spectrum, ...
     'b', 'LineWidth',2, ...
     'DisplayName','Broadband Noise (1.5f0 – 3f0)');

% Noise floor line
yline(ax, broadband_noise_avg_power_dB, 'k--', ...
      'LineWidth',2, ...
      'DisplayName',sprintf('Broadband Noise Floor = %.2f dB', broadband_noise_avg_power_dB));

% --------------------------------------------------
% 3. Fundamental + Harmonics
% --------------------------------------------------
plot(ax, f_norm(harmonic_indices), ...
     20*log10(harmonic_amplitudes), ...
     'ro', ...
     'MarkerFaceColor','r', ...
     'MarkerSize',10, ...
     'LineWidth',1.5, ...
     'DisplayName','Fundamental (f_0) and Harmonics');

% --------------------------------------------------
% 4. Ultraharmonics (Dark Green)
% --------------------------------------------------
dark_green = [0 0.5 0];

plot(ax, f_norm(ultraharmonic_indices), ...
     20*log10(ultraharmonic_amplitudes), ...
     '^', ...
     'Color', dark_green, ...
     'MarkerFaceColor', dark_green, ...
     'MarkerSize',10, ...
     'LineWidth',1.5, ...
     'DisplayName','Ultraharmonics');

% --------------------------------------------------
% 5. Annotations
% --------------------------------------------------
label_offset = 0.06 * max(spectrum_dB);
annotation_fontsize = 14;

for k = 1:n_harmonics
    if k == 1
        label = 'f_0';
    else
        label = sprintf('H%d', k-1);
    end
    
    text(f_norm(harmonic_indices(k)), ...
         20*log10(harmonic_amplitudes(k)) + label_offset, ...
         label, ...
         'Color','r', ...
         'FontWeight','bold', ...
         'FontSize',annotation_fontsize, ...
         'HorizontalAlignment','center');
end

for k = 1:length(uh_orders)
    text(f_norm(ultraharmonic_indices(k)), ...
         20*log10(ultraharmonic_amplitudes(k)) + label_offset, ...
         sprintf('U%d', k), ...
         'Color', dark_green, ...
         'FontWeight','bold', ...
         'FontSize',annotation_fontsize, ...
         'HorizontalAlignment','center');
end

% --------------------------------------------------
% Formatting: Title, Axis Labels, Legend
% --------------------------------------------------
title('Normalized Cavitation Spectrum Showing Harmonics, Ultraharmonics, and Broadband Noise (1.5f0–3f0)', ...
      'FontSize',16,'FontWeight','bold')

xlabel('Normalized Frequency  ( f / f_0 ) (MPa)', ...
       'FontSize',15,'FontWeight','bold')

ylabel('Acoustic Power Spectral Density (dB)', ...
       'FontSize',15,'FontWeight','bold')

set(ax,'FontSize',13,'LineWidth',1.5)
xlim([0 6])
grid on
box on

legend('Location','northeast','FontSize',12)
