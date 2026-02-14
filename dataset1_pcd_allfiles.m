%% in-vitro PCD Dataset — Batch Harmonic Analysis 
%% PCD Dataset — Batch Harmonic & Ultraharmonic Power Extraction + Broadband Noise
clear; clc; close all;

%% Anchor to THIS script's folder (data/)
scriptDir = fileparts(mfilename('fullpath'));
cd(scriptDir);
fprintf('Script directory: %s\n', scriptDir)

%% Add repo to path
repo_path = fullfile(scriptDir, '04_Projects', '17_FUS_Instruments_Interface');
if exist(repo_path, 'dir')
    addpath(genpath(repo_path));
    fprintf('Added repository to path\n');
else
    warning('Repo not found — continuing anyway');
end

%% Target ONLY PCD data
pcd_root = fullfile(scriptDir, 'pcd');
assert(exist(pcd_root, 'dir') == 7, 'PCD folder not found');

fprintf('PCD root: %s\n', pcd_root)

%% Analysis parameters
params.drive_frequency = 0.5e6;        % Hz
params.n_harmonics = 4;                % f0 + H1–H3
params.uh_orders = [1.5 2.5 3.5];      % ultraharmonics
params.bin_range = 1000;               % bins around peak
params.noise_exclusion_range = 2000;   % bins to exclude around peaks for noise

%% Find ALL PCD .mat files (recursive)
mat_files = dir(fullfile(pcd_root, '**', '*.mat'));

% Optional: alternate files if pressure ramps every other file
mat_files = mat_files(1:2:end);

N = length(mat_files);
fprintf('Found %d PCD files\n', N);

assert(N > 0, 'No .mat files found in pcd directory');

%% Preallocate results
Results(N) = struct( ...
    'filename', [], ...
    'folder', [], ...
    'harmonic_power', [], ...
    'harmonic_power_dB', [], ...
    'ultraharmonic_power', [], ...
    'ultraharmonic_power_dB', [], ...
    'harmonic_indices', [], ...
    'ultraharmonic_indices', [], ...
    'broadband_noise_power', [], ...
    'broadband_noise_power_dB', [], ...
    'noise_percentage', [] );

%% ================= MAIN LOOP =================
for i = 1:N
    file_path = fullfile(mat_files(i).folder, mat_files(i).name);
    fprintf('\n[%d/%d] %s\n', i, N, file_path);

    reader = load(file_path);

    % --- Welch spectrum ---
    [pxx, f] = compute_pwelch(reader.B, reader.Tinterval);

    % --- Find harmonics ---
    [harm_idx, uh_idx] = find_harmonics(f, pxx, params);

    % --- Integrate power ---
    [Hpow, Hpow_dB] = integrate_band_power(pxx, f, harm_idx, params.bin_range);
    [Upow, Upow_dB] = integrate_band_power(pxx, f, uh_idx, params.bin_range);
    
    % --- Calculate broadband noise ---
    [BBnoise, BBnoise_dB, noise_pct] = calculate_broadband_noise(pxx, f, harm_idx, uh_idx, params);

    % --- Store ---
    Results(i).filename = mat_files(i).name;
    Results(i).folder = mat_files(i).folder;
    Results(i).harmonic_power = Hpow;
    Results(i).harmonic_power_dB = Hpow_dB;
    Results(i).ultraharmonic_power = Upow;
    Results(i).ultraharmonic_power_dB = Upow_dB;
    Results(i).harmonic_indices = harm_idx;
    Results(i).ultraharmonic_indices = uh_idx;
    Results(i).broadband_noise_power = BBnoise;
    Results(i).broadband_noise_power_dB = BBnoise_dB;
    Results(i).noise_percentage = noise_pct;
    
    fprintf('  Broadband Noise: %.2f dB (%.1f%% of spectrum used)\n', BBnoise_dB, noise_pct);
end

%% Convert results to flat table (power + index + broadband noise)

% Extract arrays
H  = vertcat(Results.harmonic_power_dB);        % N × 4
U  = vertcat(Results.ultraharmonic_power_dB);   % N × 3
Hi = vertcat(Results.harmonic_indices);         % N × 4
Ui = vertcat(Results.ultraharmonic_indices);    % N × 3
BB = [Results.broadband_noise_power_dB]';       % N × 1

% Build table
ResultTable = table( ...
    {Results.filename}', ...
    Hi(:,1), Hi(:,2), Hi(:,3), Hi(:,4), ...
    H(:,1),  H(:,2),  H(:,3),  H(:,4), ...
    Ui(:,1), Ui(:,2), Ui(:,3), ...
    U(:,1),  U(:,2),  U(:,3), ...
    BB, ...
    'VariableNames', { ...
        'File', ...
        'f0_idx','H1_idx','H2_idx','H3_idx', ...
        'f0_dB','H1_dB','H2_dB','H3_dB', ...
        'U1_idx','U2_idx','U3_idx', ...
        'U1_dB','U2_dB','U3_dB', ...
        'Broadband_Noise_dB' } );

disp(ResultTable)

% Export to CSV
writetable(ResultTable, 'PCD_Harmonic_Features.csv');

% Save outputs
save('PCD_Features_AllFiles.mat', 'Results', 'ResultTable')

fprintf('\n✓ PCD feature extraction complete (with broadband noise)\n');


%% Functions
function [pxx, f] = compute_pwelch(signal, Tinterval)
% Computes Welch power spectrum using dataset sampling interval

[pxx, f] = pwelch(signal, [], [], [], 1 / Tinterval);

end

function [harmonic_indices, ultraharmonic_indices] = find_harmonics(f, pxx, params)

WelchSpec = sqrt(pxx);

% --- Harmonics ---
harmonic_freqs = params.drive_frequency * (1:params.n_harmonics);
harmonic_indices = zeros(1, params.n_harmonics);

for k = 1:params.n_harmonics
    [~, harmonic_indices(k)] = min(abs(f - harmonic_freqs(k)));
end

% --- Ultraharmonics ---
uh_freqs = params.uh_orders * params.drive_frequency;
ultraharmonic_indices = zeros(1, length(uh_freqs));

df = f(2) - f(1);
search_range = params.bin_range * df;

for k = 1:length(uh_freqs)
    idx_window = find(f >= uh_freqs(k) - search_range & ...
                      f <= uh_freqs(k) + search_range);

    if isempty(idx_window)
        ultraharmonic_indices(k) = NaN;
        continue
    end

    [~, max_idx] = max(WelchSpec(idx_window));
    ultraharmonic_indices(k) = idx_window(max_idx);
end

end

function [power_lin, power_dB] = integrate_band_power(pxx, f, indices, bin_range)
% Integrates power around frequency peaks

df = f(2) - f(1);

power_lin = zeros(1, length(indices));
power_dB  = zeros(1, length(indices));

for k = 1:length(indices)
    if isnan(indices(k))
        power_lin(k) = NaN;
        power_dB(k)  = NaN;
        continue
    end

    startBin = max(1, indices(k) - bin_range);
    endBin   = min(length(pxx), indices(k) + bin_range);

    power_lin(k) = sum(pxx(startBin:endBin)) * df;
    power_dB(k)  = 10*log10(power_lin(k));
end

end

function [noise_power, noise_power_dB, noise_percentage] = calculate_broadband_noise(pxx, f, harm_idx, uh_idx, params)
% Calculate broadband noise by excluding regions around harmonics and ultraharmonics
% Returns average power per bin in noise regions

% Create mask for all frequency bins (start with all bins included)
noise_mask = true(size(pxx));

% Combine all peak indices
all_peak_idx = [harm_idx, uh_idx];
all_peak_idx = all_peak_idx(~isnan(all_peak_idx));  % Remove NaN values

% Exclude regions around each peak
for k = 1:length(all_peak_idx)
    startBin = max(1, all_peak_idx(k) - params.noise_exclusion_range);
    endBin   = min(length(pxx), all_peak_idx(k) + params.noise_exclusion_range);
    noise_mask(startBin:endBin) = false;
end

% Get noise bins (power spectral density values)
noise_bins_psd = pxx(noise_mask);

if ~isempty(noise_bins_psd)
    % Calculate power of each noise bin
    df = f(2) - f(1);  % frequency resolution
    noise_bins_power = noise_bins_psd * df;  % convert PSD to power for each bin
    
    % Take the average of the power values
    noise_power = mean(noise_bins_power);
    noise_power_dB = 10*log10(noise_power);
    
    % Calculate percentage of spectrum used for noise
    noise_percentage = 100 * sum(noise_mask) / length(pxx);
else
    noise_power = NaN;
    noise_power_dB = NaN;
    noise_percentage = 0;
    warning('No frequency bins available for noise calculation');
end

end
