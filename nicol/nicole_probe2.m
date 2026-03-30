clear; clc; close all;
%% =========================
% PATHS
%% =========================
function_path = 'C:\Users\njcho\OneDrive\Documents\CAPSTONE\0_04_Projects\04_Projects';
probe_path    = "C:\Users\njcho\OneDrive\Documents\CAPSTONE\1_and_2_2026-01-16_single_tube_sweep\probe";
out_csv       = fullfile("C:\Users\njcho\OneDrive\Documents\CAPSTONE\1_and_2_2026-01-16_single_tube_sweep", 'probe_allfolders.csv');

addpath(function_path);

if ~exist(probe_path, 'dir')
    error('Probe folder not found: %s', probe_path);
end

%% =========================
% SETTINGS
%% =========================
f0 = 0.5e6;                % 500 kHz
probe_fs_fallback = 20e6;  % used only if fs is not available
middle_channel = 64;       % use only this channel

harmonics = 6:12;
ultras    = 6.5:1:12.5;

search_half_width_Hz = 0.08e6;
bb_low_Hz            = 9.0 * f0;
bb_high_Hz           = 10.5 * f0;
bb_exclusion_Hz      = 0.08e6;

ambientKeywords = {'ambient','control','baseline','nobubble','no_bubble','no-bubble'};

%% =========================
% FILE LIST
%% =========================
all_files = dir(fullfile(probe_path, '*.pacq*'));

if isempty(all_files)
    error('No .pacq files found in: %s', probe_path);
end

[~, order] = sort({all_files.name});
all_files = all_files(order);

keep_mask = false(numel(all_files),1);
isAmbient = false(numel(all_files),1);

for i = 1:numel(all_files)
    fname = string(all_files(i).name);
    fnameLower = lower(fname);

    % Only PCI files are allowed
    if ~startsWith(fname, "PCI_")
        continue;
    end

    % Regular experiment files:
    % PCI_..._01.pacq or PCI_..._02.pacq
    isRegular = contains(fnameLower, '_01.pacq') || contains(fnameLower, '_02.pacq');

    % Ambient PCI files only
    ambientHit = false;
    for k = 1:numel(ambientKeywords)
        if contains(fnameLower, ambientKeywords{k})
            ambientHit = true;
            break;
        end
    end

    if ambientHit
        keep_mask(i) = true;
        isAmbient(i) = true;
    elseif isRegular
        keep_mask(i) = true;
    end
end

data_files = all_files(keep_mask);
data_isAmbient = isAmbient(keep_mask);

if isempty(data_files)
    error('No matching PCI probe files found.');
end

fprintf('Found %d selected PCI probe files:\n', numel(data_files));
for i = 1:numel(data_files)
    if data_isAmbient(i)
        fprintf('  %s  [ambient]\n', data_files(i).name);
    else
        fprintf('  %s\n', data_files(i).name);
    end
end
fprintf('\n');

%% =========================
% COUNT TOTAL FRAMES
%% =========================
total_rows = 0;
for i = 1:numel(data_files)
    probe_data = read_probe_file(fullfile(data_files(i).folder, data_files(i).name), probe_fs_fallback);
    total_rows = total_rows + size(probe_data.rf_data, 3);
end

fprintf('Total output rows: %d\n\n', total_rows);

%% =========================
% PREALLOCATE
%% =========================
template = make_output_template(harmonics, ultras);
rows = repmat(template, total_rows, 1);

%% =========================
% MAIN LOOP
%% =========================
row_idx = 1;

for i = 1:numel(data_files)

    file_path = fullfile(data_files(i).folder, data_files(i).name);
    probe_data = read_probe_file(file_path, probe_fs_fallback);

    rf_data = probe_data.rf_data;
    fs      = probe_data.fs;
    nFrames = size(rf_data, 3);
    nCh     = size(rf_data, 2);

    if nCh < middle_channel
        error('File %s has only %d channels. Cannot use channel %d.', ...
            data_files(i).name, nCh, middle_channel);
    end

    if data_isAmbient(i)
        sourceType = "ambient";
    else
        sourceType = "experiment";
    end

    fprintf('Processing %s (%d frames) [%s] using channel %d\n', ...
        data_files(i).name, nFrames, sourceType, middle_channel);

    for frame_idx = 1:nFrames

        % Use ONLY channel 64
        sig = double(rf_data(:, middle_channel, frame_idx));
        sig = sig(:);

        [pxx_use, f_ref] = pwelch(sig, [], [], [], fs);
        pxx_use = pxx_use(:);
        f_ref   = f_ref(:);

        feat = template;
        feat.SourceType = sourceType;
        feat.filename = string(data_files(i).name);
        feat.frame = frame_idx;

        peak_freqs = [];

        % Harmonics
        for h = harmonics
            target_f = h * f0;
            [pk_freq, pk_db, pk_power] = extract_peak_band(f_ref, pxx_use, target_f, search_half_width_Hz);

            feat.(sprintf('H%d_freq', h))  = pk_freq;
            feat.(sprintf('H%d_db', h))    = pk_db;
            feat.(sprintf('H%d_power', h)) = pk_power;

            if ~isnan(pk_freq)
                peak_freqs(end+1,1) = pk_freq; %#ok<AGROW>
            end
        end

        % Ultraharmonics
        for u = ultras
            target_f = u * f0;
            u_label = ultra_label(u);   % valid MATLAB name, e.g. U6p5

            [pk_freq, pk_db, pk_power] = extract_peak_band(f_ref, pxx_use, target_f, search_half_width_Hz);

            feat.(sprintf('%s_freq',  u_label)) = pk_freq;
            feat.(sprintf('%s_db',    u_label)) = pk_db;
            feat.(sprintf('%s_power', u_label)) = pk_power;

            if ~isnan(pk_freq)
                peak_freqs(end+1,1) = pk_freq; %#ok<AGROW>
            end
        end

        % Broadband noise
        [bb_db, bb_power] = extract_broadband_noise( ...
            f_ref, pxx_use, bb_low_Hz, bb_high_Hz, peak_freqs, bb_exclusion_Hz);

        feat.broadband_noise_db    = bb_db;
        feat.broadband_noise_power = bb_power;

        rows(row_idx) = feat;
        row_idx = row_idx + 1;
    end
end

rows = rows(1:row_idx-1);

%% =========================
% SAVE TABLE
%% =========================
T = struct2table(rows);
T = rename_ultra_columns_for_export(T, ultras);

disp(T(1:min(5,height(T)), :));
writetable(T, out_csv);

fprintf('\nSaved output to:\n%s\n', out_csv);

%% ========================================================================
% LOCAL FUNCTIONS
%% ========================================================================

function probe_data = read_probe_file(file_path, fs_fallback)

    acq = imaging.common.acquisition_reader(string(file_path));

    N_t      = acq.rfdata_dimensions(1);
    N_ch     = acq.rfdata_dimensions(2);
    N_frames = acq.number_of_frames;

    rf_data = zeros(N_t, N_ch, N_frames, 'double');
    for k = 1:N_frames
        rf_data(:,:,k) = double(acq.data{k});
    end

    fs = [];
    try
        if ~isempty(acq.data_params)
            if isstruct(acq.data_params) && isfield(acq.data_params, 'fs')
                fs = acq.data_params.fs;
            elseif isobject(acq.data_params)
                fs = acq.data_params.fs;
            end
        end
    catch
        fs = [];
    end

    if isempty(fs)
        fs = fs_fallback;
    end

    probe_data = struct();
    probe_data.rf_data = rf_data;
    probe_data.fs = fs;
end

function feat = make_output_template(harmonics, ultras)

    feat = struct();

    feat.SourceType = "";

    for h = harmonics
        feat.(sprintf('H%d_freq', h))  = NaN;
        feat.(sprintf('H%d_db', h))    = NaN;
        feat.(sprintf('H%d_power', h)) = NaN;
    end

    for u = ultras
        u_label = ultra_label(u);
        feat.(sprintf('%s_freq',  u_label)) = NaN;
        feat.(sprintf('%s_db',    u_label)) = NaN;
        feat.(sprintf('%s_power', u_label)) = NaN;
    end

    feat.broadband_noise_db    = NaN;
    feat.broadband_noise_power = NaN;
    feat.filename = "";
    feat.frame    = NaN;
end

function lbl = ultra_label(u)
    lbl = ['U' strrep(sprintf('%.1f', u), '.', 'p')];
end

function T = rename_ultra_columns_for_export(T, ultras)
    names = T.Properties.VariableNames;

    for u = ultras
        internal = ultra_label(u);      % U6p5
        external = sprintf('U%.1f', u); % U6.5

        names = strrep(names, [internal '_freq'],  [external '_freq']);
        names = strrep(names, [internal '_db'],    [external '_db']);
        names = strrep(names, [internal '_power'], [external '_power']);
    end

    T.Properties.VariableNames = names;
end

function [peak_freq, peak_db, band_power] = extract_peak_band(f, pxx, target_f, half_width)

    f = f(:);
    pxx = pxx(:);
    pxx_dB = 10 * log10(pxx + eps);

    idx = find(f >= (target_f - half_width) & f <= (target_f + half_width));

    if isempty(idx)
        peak_freq = NaN;
        peak_db = NaN;
        band_power = NaN;
        return;
    end

    [peak_db, rel_idx] = max(pxx_dB(idx));
    peak_idx = idx(rel_idx);
    peak_freq = f(peak_idx);

    band_power = safe_band_power(f(idx), pxx(idx));
end

function [bb_db, bb_power] = extract_broadband_noise(f, pxx, low_Hz, high_Hz, peak_freqs, exclusion_Hz)

    f = f(:);
    pxx = pxx(:);

    mask = (f >= low_Hz) & (f <= high_Hz);

    for k = 1:numel(peak_freqs)
        mask = mask & (abs(f - peak_freqs(k)) > exclusion_Hz);
    end

    f_noise = f(mask);
    pxx_noise = pxx(mask);

    if isempty(pxx_noise)
        bb_db = NaN;
        bb_power = NaN;
        return;
    end

    bb_db = 10 * log10(mean(pxx_noise) + eps);
    bb_power = safe_band_power(f_noise, pxx_noise);
end

function val = safe_band_power(f_segment, pxx_segment)

    f_segment = f_segment(:);
    pxx_segment = pxx_segment(:);

    if isempty(f_segment) || isempty(pxx_segment)
        val = NaN;
        return;
    end

    if numel(f_segment) ~= numel(pxx_segment)
        val = NaN;
        return;
    end

    if numel(f_segment) == 1
        val = pxx_segment(1);
        return;
    end

    val = trapz(f_segment, pxx_segment);
end