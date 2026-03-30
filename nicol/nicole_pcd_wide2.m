clear; clc;
%% ---------------- USER SETTINGS ----------------
% Use parent pcd folder so all sibling folders are included
pcdRoot = 'C:\Users\njcho\OneDrive\Documents\CAPSTONE\1_and_2_2026-01-16_single_tube_sweep\pcd';

saveCSV = fullfile('C:\Users\njcho\OneDrive\Documents\CAPSTONE\1_and_2_2026-01-16_single_tube_sweep', 'pcd_allfolders.csv');

f0 = 0.5e6;                % fundamental frequency in Hz
searchWidthHz = 0.10e6;    % UH peak search half-width
bbWindowHz    = 0.05e6;    % broadband averaging half-width
guardWidthHz  = 0.005e6;   % exclude around UH peak

harmonicMults = [1 2 3 4];     % f0, 2f0, 3f0, 4f0
uhMults       = [1.5 2.5 3.5]; % ultraharmonics

% Ambient inclusion
includeAmbient = true;

% Folder-name keywords that count as ambient/control folders
ambientKeywords = {'ambient','control','baseline','nobubble','no_bubble','no-bubble'};

%% ---------------- FIND MATCHING SUBFOLDERS ----------------
allEntries = dir(pcdRoot);
isSub = [allEntries.isdir];
subFolders = allEntries(isSub);
subFolders = subFolders(~ismember({subFolders.name}, {'.','..'}));

folderNames = {subFolders.name};

keepFolder = false(size(folderNames));
isAmbientFolder = false(size(folderNames));

for i = 1:numel(folderNames)
    thisName = folderNames{i};

    % Regular experiment folders like:
    % 1m_x_01, 1m_x_02, 100x_01, 100x_02, 1000x_01, 1000x_02
    isRegular = ~isempty(regexp(thisName, '^\d+m?_x?_\d+$|^\d+x_\d+$', 'once'));

    % Ambient folders by keyword
    ambientHit = false;
    for k = 1:numel(ambientKeywords)
        if contains(lower(thisName), lower(ambientKeywords{k}))
            ambientHit = true;
            break;
        end
    end

    if isRegular
        keepFolder(i) = true;
    elseif includeAmbient && ambientHit
        keepFolder(i) = true;
        isAmbientFolder(i) = true;
    end
end

subFolders = subFolders(keepFolder);
isAmbientFolder = isAmbientFolder(keepFolder);

if isempty(subFolders)
    error('No matching PCD subfolders found inside: %s', pcdRoot);
end

[~, ord] = sort({subFolders.name});
subFolders = subFolders(ord);

% Reorder ambient flags consistently
tmpNames = {subFolders.name};
ambientFlagsSorted = false(size(tmpNames));
for i = 1:numel(tmpNames)
    for j = 1:numel(folderNames)
        if strcmp(tmpNames{i}, folderNames{j})
            ambientFlagsSorted(i) = ismember(lower(folderNames{j}), lower(folderNames(j))) && false;
        end
    end
end

% Recompute ambient flags directly from sorted names
ambientFlagsSorted = false(numel(subFolders),1);
for i = 1:numel(subFolders)
    thisName = lower(subFolders(i).name);
    for k = 1:numel(ambientKeywords)
        if contains(thisName, lower(ambientKeywords{k}))
            ambientFlagsSorted(i) = true;
            break;
        end
    end
end

fprintf('Folders to process:\n');
for i = 1:numel(subFolders)
    if ambientFlagsSorted(i)
        fprintf('  %s  [ambient]\n', subFolders(i).name);
    else
        fprintf('  %s\n', subFolders(i).name);
    end
end
fprintf('\n');

%% ---------------- GATHER EVEN-NUMBERED MAT FILES ----------------
fileList = struct('folder', {}, 'name', {}, 'fullpath', {}, 'group', {}, 'isAmbient', {});

for i = 1:numel(subFolders)
    thisFolder = fullfile(subFolders(i).folder, subFolders(i).name);
    mats = dir(fullfile(thisFolder, '*.mat'));

    if isempty(mats)
        fprintf('No .mat files in %s\n', subFolders(i).name);
        continue;
    end

    [~, mOrd] = sort({mats.name});
    mats = mats(mOrd);

    for j = 1:numel(mats)
        tok = regexp(mats(j).name, '_(\d+)\.mat$', 'tokens');

        % Skip files that do not end with _NN.mat
        if isempty(tok)
            continue;
        end

        fileNum = str2double(tok{1}{1});

        % only even-numbered filenames
        if mod(fileNum, 2) ~= 0
            continue;
        end

        fileList(end+1).folder    = mats(j).folder; %#ok<SAGROW>
        fileList(end).name        = mats(j).name;
        fileList(end).fullpath    = fullfile(mats(j).folder, mats(j).name);
        fileList(end).group       = subFolders(i).name;
        fileList(end).isAmbient   = ambientFlagsSorted(i);
    end
end

if isempty(fileList)
    error('No even-numbered .mat files found in the selected folders.');
end

fprintf('Total even-numbered .mat files selected: %d\n\n', numel(fileList));

%% ---------------- PREALLOCATE ----------------
nFiles = numel(fileList);

Folder     = strings(nFiles,1);
File       = strings(nFiles,1);
SourceType = strings(nFiles,1);

f0_idx = zeros(nFiles,1);
H1_idx = zeros(nFiles,1);
H2_idx = zeros(nFiles,1);
H3_idx = zeros(nFiles,1);

f0_dB = zeros(nFiles,1);
H1_dB = zeros(nFiles,1);
H2_dB = zeros(nFiles,1);
H3_dB = zeros(nFiles,1);

U1_idx = zeros(nFiles,1);
U2_idx = zeros(nFiles,1);
U3_idx = zeros(nFiles,1);

U1_dB = zeros(nFiles,1);
U2_dB = zeros(nFiles,1);
U3_dB = zeros(nFiles,1);

U1_BB_dB = zeros(nFiles,1);
U2_BB_dB = zeros(nFiles,1);
U3_BB_dB = zeros(nFiles,1);

%% ---------------- MAIN LOOP ----------------
for k = 1:nFiles
    S = load(fileList(k).fullpath);

    if ~isfield(S, 'B')
        error('File %s does not contain channel B.', fileList(k).name);
    end
    if ~isfield(S, 'Tinterval')
        error('File %s does not contain Tinterval.', fileList(k).name);
    end

    x  = double(S.B(:));
    Fs = 1 / double(S.Tinterval);

    % PWELCH (replaces FFT block)
    window   = hamming(1024);
    noverlap = 512;
    nfft     = 2048;

    [pxx, f] = pwelch(x, window, noverlap, nfft, Fs);
    spec_dB = 10 * log10(pxx + eps);

    Folder(k) = string(fileList(k).group);
    File(k)   = string(fileList(k).name);

    if fileList(k).isAmbient
        SourceType(k) = "ambient";
    else
        SourceType(k) = "experiment";
    end

    %% Harmonics
    harmIdx = zeros(1, numel(harmonicMults));
    harmVal = zeros(1, numel(harmonicMults));

    for h = 1:numel(harmonicMults)
        targetFreq = harmonicMults(h) * f0;
        [harmIdx(h), harmVal(h)] = nearestBinValue(f, spec_dB, targetFreq);
    end

    f0_idx(k) = harmIdx(1);
    H1_idx(k) = harmIdx(2);
    H2_idx(k) = harmIdx(3);
    H3_idx(k) = harmIdx(4);

    f0_dB(k) = harmVal(1);
    H1_dB(k) = harmVal(2);
    H2_dB(k) = harmVal(3);
    H3_dB(k) = harmVal(4);

    %% Ultraharmonics
    uhIdx = zeros(1, numel(uhMults));
    uhVal = zeros(1, numel(uhMults));
    uhBB  = zeros(1, numel(uhMults));

    for u = 1:numel(uhMults)
        targetFreq = uhMults(u) * f0;

        [uhIdx(u), uhVal(u)] = localPeakInWindow(f, spec_dB, targetFreq, searchWidthHz);

        peakFreq = f(uhIdx(u));
        uhBB(u) = broadbandAroundPeak(f, spec_dB, peakFreq, bbWindowHz, guardWidthHz);
    end

    U1_idx(k) = uhIdx(1);
    U2_idx(k) = uhIdx(2);
    U3_idx(k) = uhIdx(3);

    U1_dB(k) = uhVal(1);
    U2_dB(k) = uhVal(2);
    U3_dB(k) = uhVal(3);

    U1_BB_dB(k) = uhBB(1);
    U2_BB_dB(k) = uhBB(2);
    U3_BB_dB(k) = uhBB(3);

    fprintf('Processed %d/%d: %s | %s | %s\n', ...
        k, nFiles, SourceType(k), fileList(k).group, fileList(k).name);
end

%% ---------------- BUILD TABLE ----------------
T = table( ...
    SourceType, Folder, File, ...
    f0_idx, H1_idx, H2_idx, H3_idx, ...
    f0_dB, H1_dB, H2_dB, H3_dB, ...
    U1_idx, U2_idx, U3_idx, ...
    U1_dB, U2_dB, U3_dB, ...
    U1_BB_dB, U2_BB_dB, U3_BB_dB);

disp(T(1:min(10,height(T)), :));
writetable(T, saveCSV);

fprintf('\nFinished.\nSaved to:\n%s\n', saveCSV);

%% ---------------- LOCAL FUNCTIONS ----------------
function [idx, val_dB] = nearestBinValue(f, spec_dB, targetFreq)
    [~, idx] = min(abs(f - targetFreq));
    val_dB = spec_dB(idx);
end

function [idxPeak, peakdB] = localPeakInWindow(f, spec_dB, targetFreq, halfWidthHz)
    mask = (f >= targetFreq - halfWidthHz) & (f <= targetFreq + halfWidthHz);
    idxLocal = find(mask);

    if isempty(idxLocal)
        [~, idxPeak] = min(abs(f - targetFreq));
        peakdB = spec_dB(idxPeak);
        return
    end

    [~, relMax] = max(spec_dB(idxLocal));
    idxPeak = idxLocal(relMax);
    peakdB = spec_dB(idxPeak);
end

function bb_dB = broadbandAroundPeak(f, spec_dB, peakFreq, bbHalfWidthHz, guardHalfWidthHz)
    outerMask = (f >= peakFreq - bbHalfWidthHz) & (f <= peakFreq + bbHalfWidthHz);
    guardMask = (f >= peakFreq - guardHalfWidthHz) & (f <= peakFreq + guardHalfWidthHz);

    bbMask = outerMask & ~guardMask;
    vals = spec_dB(bbMask);

    if isempty(vals)
        bb_dB = NaN;
    else
        bb_dB = mean(vals, 'omitnan');
    end
end