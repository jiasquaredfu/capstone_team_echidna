clear; clc;
%% ---------------- PATHS ----------------
pcd_csv   = 'C:\Users\njcho\OneDrive\Documents\CAPSTONE\1_and_2_2026-01-16_single_tube_sweep\pcd_allfolders.csv';
probe_csv = 'C:\Users\njcho\OneDrive\Documents\CAPSTONE\1_and_2_2026-01-16_single_tube_sweep\probe_allfolders.csv';

out_csv   = 'C:\Users\njcho\OneDrive\Documents\CAPSTONE\1_and_2_2026-01-16_single_tube_sweep\matched_pcd_probe.csv';

%% ---------------- READ TABLES ----------------
PCD   = readtable(pcd_csv);
PROBE = readtable(probe_csv);

%% ---------------- STANDARDIZE STRING COLUMNS ----------------
if ismember('Folder', PCD.Properties.VariableNames)
    PCD.Folder = string(PCD.Folder);
end
if ismember('File', PCD.Properties.VariableNames)
    PCD.File = string(PCD.File);
end
if ismember('filename', PROBE.Properties.VariableNames)
    PROBE.filename = string(PROBE.filename);
end

%% ---------------- BUILD MATCH KEYS ----------------
% PCD condition from folder name
PCD.condition = string(PCD.Folder);

% PCD frame from even-numbered .mat filename
% Example: ..._02.mat -> frame 1, ..._04.mat -> frame 2
PCD.frame = nan(height(PCD),1);

for i = 1:height(PCD)
    fname = PCD.File(i);
    tok = regexp(fname, '_(\d+)\.mat$', 'tokens', 'once');

    if ~isempty(tok)
        fileNum = str2double(tok{1});
        PCD.frame(i) = fileNum / 2;
    end
end

% Probe condition from filename
% Examples:
% BMode_1000x_01.pacq -> 1000x_01
% PCI_1000x_01.pacq   -> 1000x_01
PROBE.condition = strings(height(PROBE),1);

for i = 1:height(PROBE)
    fname = PROBE.filename(i);
    tok = regexp(fname, '(?:BMode_|PCI_)(.+)\.pacq$', 'tokens', 'once');

    if ~isempty(tok)
        PROBE.condition(i) = string(tok{1});
    end
end

% Probe frame should already exist
PROBE.frame = double(PROBE.frame);

%% ---------------- CLEAN BAD ROWS ----------------
PCD = PCD(~ismissing(PCD.condition) & PCD.condition ~= "" & ~isnan(PCD.frame), :);
PROBE = PROBE(~ismissing(PROBE.condition) & PROBE.condition ~= "" & ~isnan(PROBE.frame), :);

%% ---------------- PREFIX COLUMN NAMES ----------------
pcdNames = PCD.Properties.VariableNames;
probeNames = PROBE.Properties.VariableNames;

PCD.Properties.VariableNames = strcat("pcd_", pcdNames);
PROBE.Properties.VariableNames = strcat("probe_", probeNames);

%% ---------------- JOIN ----------------
Matched = innerjoin( ...
    PCD, PROBE, ...
    'LeftKeys',  {'pcd_condition','pcd_frame'}, ...
    'RightKeys', {'probe_condition','probe_frame'} );
%% ---------------- OPTIONAL SORT ----------------
Matched = sortrows(Matched);

%% ---------------- SAVE ----------------
writetable(Matched, out_csv);

fprintf('Matched rows: %d\n', height(Matched));
fprintf('Saved matched table to:\n%s\n\n', out_csv);

disp(Matched(1:min(10,height(Matched)), :));