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

%raw signal just to see
fig = figure('Name','Raw Signals','NumberTitle','off','Position',[100 100 1000 600]);
tgroup = uitabgroup(fig);

% --- Tab 1: A signal ---
tab1 = uitab(tgroup,'Title','A signal');
ax1 = axes('Parent', tab1);
plot(ax1, t, reader.A, 'b');
title(ax1,'Raw A Signal')
xlabel(ax1,'Time (s)')
ylabel(ax1,'Amplitude')
grid(ax1,'on')

% --- Tab 2: B signal ---
tab2 = uitab(tgroup,'Title','B signal');
ax2 = axes('Parent', tab2);
plot(ax2, t, reader.B, 'r');
title(ax2,'Raw B Signal')
xlabel(ax2,'Time (s)')
ylabel(ax2,'Amplitude')
grid(ax2,'on')


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

% number of points and sampling frequency
N = length(reader.B);
Fs = 1 / reader.Tinterval;

Y = fft(reader.B);
% two-sided spectrum magnitude
P2 = abs(Y / N);
%single-sides specturm magnitude
P1 = P2(1:floor(N/2)+1);
%frequency vector corresponding to P1
f = Fs * (0:floor(N/2)) / N;
%0.5 MHz fundamental
drive_frequency = 0.5e6; 
%number of harmonics we are looking for
n_harmonics = 3;            
%array of harmonics
harmonic_freqs = drive_frequency * (1:n_harmonics); 
harmonic_indices = zeros(1, n_harmonics);
for k = 1:n_harmonics
    [~, idx] = min(abs(f - harmonic_freqs(k)));  
    harmonic_indices(k) = idx;
end
harmonic_amplitudes = P1(harmonic_indices);

% --- Ultraharmonics ---
uh_multiples = 3:2:(2*n_harmonics+1); 
ultraharmonic_freqs = uh_multiples * drive_frequency / 2;

ultraharmonic_indices = zeros(1, length(ultraharmonic_freqs));
ultraharmonic_amplitudes = zeros(1, length(ultraharmonic_freqs));

for k = 1:length(ultraharmonic_freqs)
    [~, idx] = min(abs(f - ultraharmonic_freqs(k)));   
    ultraharmonic_indices(k) = idx;
    ultraharmonic_amplitudes(k) = P1(idx);         
end

%labeled fft graph
tab2 = uitab(tgroup,'Title','Harmonics labeled');
ax2 = axes('Parent', tab2);
plot(f / 1e6, 20*log10(P1));
hold on
plot(f(harmonic_indices) / 1e6, 20*log10(harmonic_amplitudes), 'ro', 'MarkerFaceColor','r');
title(ax2,'FFT B signal with Harmonics')
xlabel(ax2,'Frequency (MHz)')
ylabel(ax2,' Magnitude (dB)')
xlim([0 3])
grid(ax1,'on')

plot(f(ultraharmonic_indices)/1e6, 20*log10(ultraharmonic_amplitudes), 'go', 'MarkerFaceColor','g');
label_offset = 0.05 * max(20*log10(P1));

% Harmonics
for k = 1:n_harmonics
    text(f(harmonic_indices(k))/1e6, 20*log10(harmonic_amplitudes(k)) + label_offset, ...
        sprintf('H%d', k), 'Color', 'r', 'FontWeight','bold', 'HorizontalAlignment','center');
end

% Ultraharmonics
for k = 1:length(ultraharmonic_indices)
    text(f(ultraharmonic_indices(k))/1e6, 20*log10(ultraharmonic_amplitudes(k)) + label_offset, ...
        sprintf('UH%.1f', uh_multiples(k)/2), 'Color', 'g', 'FontWeight','bold', 'HorizontalAlignment','center');
end

% --- Harmonics ---
harmonics_table = [harmonic_freqs(:), harmonic_indices(:), harmonic_amplitudes(:)];
fprintf('\n  Harmonic Frequency (Hz)   FFT Index   Amplitude\n');
for k = 1:size(harmonics_table,1)
    fprintf('  %20.2f %15d %12.5f\n', harmonics_table(k,1), harmonics_table(k,2), harmonics_table(k,3));
end

% --- Ultraharmonics ---
ultraharmonics_table = [ultraharmonic_freqs(:), ultraharmonic_indices(:), ultraharmonic_amplitudes(:)];
fprintf('\n  Ultraharmonic Frequency (Hz)   FFT Index   Amplitude\n');
for k = 1:size(ultraharmonics_table,1)
    fprintf('  %20.2f %15d %12.5f\n', ultraharmonics_table(k,1), ultraharmonics_table(k,2), ultraharmonics_table(k,3));
end

%% 
