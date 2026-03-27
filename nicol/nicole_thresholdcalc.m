clear; clc;
%% ---------------- PATHS ----------------
basePath = 'C:\Users\njcho\OneDrive\Documents\CAPSTONE\1_and_2_2026-01-16_single_tube_sweep';

in_csv  = fullfile(basePath, 'matched_pcd_probe.csv');
out_csv = fullfile(basePath, 'matched_pcd_probe.csv');

%% ---------------- READ TABLE ----------------
T = readtable(in_csv);

%% ---------------- BASIC CHECKS ----------------
requiredVars = { ...
    'pcd_SourceType', ...
    'pcd_condition', ...
    'pcd_U1_BB_dB', ...
    'pcd_U2_BB_dB', ...
    'pcd_U3_BB_dB'};

for i = 1:numel(requiredVars)
    if ~ismember(requiredVars{i}, T.Properties.VariableNames)
        error('Missing required column: %s', requiredVars{i});
    end
end

T.pcd_SourceType = lower(string(T.pcd_SourceType));
T.pcd_condition  = string(T.pcd_condition);

%% ---------------- DEFINE ONE BB VALUE PER ROW ----------------
bbMat = [T.pcd_U1_BB_dB, T.pcd_U2_BB_dB, T.pcd_U3_BB_dB];
T.row_BB_dB = max(bbMat, [], 2, 'omitnan');

%% ---------------- PREALLOCATE OUTPUT ----------------
T.bb_threshold_k1 = nan(height(T),1);
T.bb_threshold_k2 = nan(height(T),1);
T.bb_threshold_k3 = nan(height(T),1);

T.unsafe_k1 = nan(height(T),1);
T.unsafe_k2 = nan(height(T),1);
T.unsafe_k3 = nan(height(T),1);

%% ---------------- EXTRACT CONCENTRATIONS ----------------
allConditions = unique(T.pcd_condition(T.pcd_SourceType ~= "ambient"));

expConc = strings(numel(allConditions),1);
for i = 1:numel(allConditions)
    expConc(i) = extract_concentration(allConditions(i));
end

uniqueConc = unique(expConc(expConc ~= ""));

%% ---------------- AMBIENT SUBSET ----------------
A = T(T.pcd_SourceType == "ambient", :);

if isempty(A)
    error('No ambient rows found.');
end

A.ambient_conc = strings(height(A),1);
for i = 1:height(A)
    A.ambient_conc(i) = extract_concentration_from_ambient_name(A.pcd_condition(i));
end

%% ---------------- COMPUTE THRESHOLDS ----------------
for c = 1:numel(uniqueConc)
    conc = uniqueConc(c);

    % Experimental rows
    expMask = (T.pcd_SourceType ~= "ambient") & startsWith(T.pcd_condition, conc);

    % Ambient rows for this concentration (if labeled)
    ambMask_specific = A.ambient_conc == conc;

    if any(ambMask_specific)
        bbAmbient = A.row_BB_dB(ambMask_specific);
    else
        % fallback: use all ambient
        bbAmbient = A.row_BB_dB;
    end

    bbAmbient = bbAmbient(~isnan(bbAmbient));

    if isempty(bbAmbient)
        warning('No ambient data for %s', conc);
        continue;
    end

    mu  = mean(bbAmbient, 'omitnan');
    sig = std(bbAmbient, 'omitnan');

    thr1 = mu + 1*sig;
    thr2 = mu + 2*sig;
    thr3 = mu + 3*sig;

    T.bb_threshold_k1(expMask) = thr1;
    T.bb_threshold_k2(expMask) = thr2;
    T.bb_threshold_k3(expMask) = thr3;

    T.unsafe_k1(expMask) = double(T.row_BB_dB(expMask) >= thr1);
    T.unsafe_k2(expMask) = double(T.row_BB_dB(expMask) >= thr2);
    T.unsafe_k3(expMask) = double(T.row_BB_dB(expMask) >= thr3);

    fprintf('%s → mean=%.3f std=%.3f | thr1=%.3f thr2=%.3f thr3=%.3f\n', ...
        conc, mu, sig, thr1, thr2, thr3);
end

%% ---------------- REMOVE AMBIENT ROWS ----------------
T = T(T.pcd_SourceType ~= "ambient", :);

%% ---------------- REMOVE pcd_SourceType COLUMN ----------------
if ismember('pcd_SourceType', T.Properties.VariableNames)
    T.pcd_SourceType = [];
end

%% ---------------- MOVE NEW COLUMNS TO END ----------------
newCols = { ...
    'row_BB_dB', ...
    'bb_threshold_k1', 'bb_threshold_k2', 'bb_threshold_k3', ...
    'unsafe_k1', 'unsafe_k2', 'unsafe_k3'};

existingNames = T.Properties.VariableNames;
keepCols = existingNames(~ismember(existingNames, newCols));

T = T(:, [keepCols, newCols]);
%% ---------------- REMOVE UNWANTED COLUMNS ----------------
vars = string(T.Properties.VariableNames);

removeMask = false(size(vars));

for i = 1:numel(vars)
    v = vars(i);

    if v == "pcd_Folder" || ...
       v == "pcd_File" || ...
       v == "pcd_frame" || ...
       v == "probe_condition" || ...
       v == "probe_frame" || ...
       v == "probe_filename" || ...
       v == "pcd_f0_dB" || ...
       contains(v, "f0") || ...
       contains(v, "BB") || ...
       contains(v, "broadband") || ...
       contains(v, "threshold")

        removeMask(i) = true;
    end
end

T(:, removeMask) = [];
%% ---------------- SAVE ----------------
writetable(T, out_csv);

fprintf('\nSaved to:\n%s\n', out_csv);
fprintf('Final row count: %d\n', height(T));

disp(T(1:min(10,height(T)), :));

%% ---------------- LOCAL FUNCTIONS ----------------
function conc = extract_concentration(condStr)
    condStr = string(condStr);
    tok = regexp(condStr, '^([0-9]+m|[0-9]+x)', 'tokens', 'once');

    if isempty(tok)
        conc = "";
    else
        conc = string(tok{1});
    end
end

function conc = extract_concentration_from_ambient_name(condStr)
    condStr = lower(string(condStr));
    tok = regexp(condStr, '([0-9]+m|[0-9]+x)', 'tokens', 'once');

    if isempty(tok)
        conc = "";
    else
        conc = string(tok{1});
    end
end