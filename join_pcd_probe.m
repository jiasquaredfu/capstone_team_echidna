%% Join PCD and Probe feature tables
% Matches rows on: concentration + pressure_MPa
%
% PCD key  : concentration parsed from 'File' column (e.g. '1000x_01_24.mat' → '1000x_01')
%            pressure taken directly from 'Pressure_MP_10' column (already in CSV)
% Probe key: 'concentration' and 'pressure_MPa' columns written by probe_feature_extraction.m
%
% Run AFTER both CSVs exist:
%   PCD_Harmonic_Features.csv
%   Probe_Features.csv

clear; clc;

%% --- Load tables ---
pcd_file   = 'PCD_Harmonic_Features.csv';
probe_file = 'Probe_Features.csv';

assert(isfile(pcd_file),   'PCD CSV not found: %s',   pcd_file);
assert(isfile(probe_file), 'Probe CSV not found: %s', probe_file);

T_pcd   = readtable(pcd_file,   'TextType', 'string');
T_probe = readtable(probe_file, 'TextType', 'string');

fprintf('PCD rows:   %d\n', height(T_pcd));
fprintf('Probe rows: %d\n', height(T_probe));

fprintf('\nPCD columns:\n');   disp(T_pcd.Properties.VariableNames');
fprintf('\nProbe columns:\n'); disp(T_probe.Properties.VariableNames');

%% --- FIX 1: Standardize PCD pressure column name ---
if ismember('Pressure_MPa', T_pcd.Properties.VariableNames)
    T_pcd = renamevars(T_pcd, 'Pressure_MPa', 'pressure_MPa');
end

%% --- Build concentration key for PCD ---
% Example: '1000x_01_24.mat' → '1000x_01'
n_pcd    = height(T_pcd);
pcd_conc = strings(n_pcd, 1);

for i = 1:n_pcd
    [~, stem, ~] = fileparts(T_pcd.File(i));   % '1000x_01_24'
    parts        = strsplit(stem, '_');
    pcd_conc(i)  = strjoin(parts(1:end-1), '_');  % remove last token
end

T_pcd.concentration = pcd_conc;

%% --- Normalize keys in BOTH tables (CRITICAL) ---
T_pcd.concentration   = string(strtrim(T_pcd.concentration));
T_probe.concentration = string(strtrim(T_probe.concentration));

T_pcd.pressure_MPa   = round(double(T_pcd.pressure_MPa), 4);
T_probe.pressure_MPa = round(double(T_probe.pressure_MPa), 4);

%% --- Diagnostics BEFORE join ---
fprintf('\nPCD  pressure range:  %.4f - %.4f MPa\n', ...
    min(T_pcd.pressure_MPa), max(T_pcd.pressure_MPa));

fprintf('Probe pressure range: %.4f - %.4f MPa\n', ...
    min(T_probe.pressure_MPa), max(T_probe.pressure_MPa));

fprintf('\nPCD concentrations:\n');
disp(unique(T_pcd.concentration));

fprintf('Probe concentrations:\n');
disp(unique(T_probe.concentration));

fprintf('\nSample PCD keys:\n');
pcd_keys_sample = unique(T_pcd(:, {'concentration','pressure_MPa'}), 'rows');
disp(pcd_keys_sample(1:min(10,height(pcd_keys_sample)), :));

fprintf('\nSample Probe keys:\n');
probe_keys_sample = unique(T_probe(:, {'concentration','pressure_MPa'}), 'rows');
disp(probe_keys_sample(1:min(10,height(probe_keys_sample)), :));

%% --- Prefix non-key columns to avoid collisions ---
shared_keys = {'concentration', 'pressure_MPa'};

pcd_cols   = T_pcd.Properties.VariableNames;
probe_cols = T_probe.Properties.VariableNames;

for k = 1:numel(pcd_cols)
    col = pcd_cols{k};
    if ~ismember(col, shared_keys)
        T_pcd = renamevars(T_pcd, col, ['pcd_' col]);
    end
end

for k = 1:numel(probe_cols)
    col = probe_cols{k};
    if ~ismember(col, shared_keys)
        T_probe = renamevars(T_probe, col, ['probe_' col]);
    end
end

%% --- Inner join ---
T_joined = innerjoin(T_pcd, T_probe, 'Keys', shared_keys);

fprintf('\nJoined table: %d rows x %d cols\n', height(T_joined), width(T_joined));

%% --- If join fails, show why ---
if height(T_joined) == 0
    warning('No matching rows found.');

    pcd_keys   = unique(T_pcd(:,   {'concentration','pressure_MPa'}), 'rows');
    probe_keys = unique(T_probe(:, {'concentration','pressure_MPa'}), 'rows');

    fprintf('\nPCD unique keys (first 10):\n');
    disp(pcd_keys(1:min(10,height(pcd_keys)), :));

    fprintf('\nProbe unique keys (first 10):\n');
    disp(probe_keys(1:min(10,height(probe_keys)), :));

    fprintf('\nConcentration mismatches:\n');
    disp(setdiff(unique(T_pcd.concentration), unique(T_probe.concentration)));
    disp(setdiff(unique(T_probe.concentration), unique(T_pcd.concentration)));

else
    %% --- Save result ---
    outfile = 'Joined_PCD_Probe_Features.csv';
    writetable(T_joined, outfile);
    fprintf('Saved: %s\n', outfile);

    fprintf('\nFirst 3 rows:\n');
    disp(T_joined(1:min(3,height(T_joined)), :));
end