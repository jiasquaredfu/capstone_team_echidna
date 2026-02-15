clear; clc; close all;
rng(42);

% load probe features 
T = readtable('dataset2_probe_data_all_frames.csv');

fprintf('Loaded %d rows from probe_data_all_frames.csv\n', height(T));
disp(T.Properties.VariableNames');

%% build a reduced feature matrix 
%  Full extraction has H6-H20 + U6-U20 + BB = 31 features
%  too many for ~224 samples (overfitting risk)
%
%  for a reduced set: pick 4 harmonics + 4 ultraharmonics spread across
%  the L12-5 probe's sensitive bandwidth (5-12 MHz), plus BB noise
%  this selection is guided by SHAP ranking from full-feature run:
%    H8 (4.0 MHz), H12 (6.0 MHz), H16 (8.0 MHz), H18 (9.0 MHz)
%    U8 (4.25 MHz), U12 (6.25 MHz), U16 (8.25 MHz), U18 (9.25 MHz)

selected_harmonics = [8, 12, 16, 18];
selected_ultraharmonics = [8, 12, 16, 18];

col_names = {};
feature_labels = {};

% harmonic dB columns
for h = selected_harmonics
    col = sprintf('H%d_db', h);
    if ismember(col, T.Properties.VariableNames)
        col_names{end+1} = col;
        feature_labels{end+1} = sprintf('H%d (%.1f MHz)', h, h*0.5);
    else
        fprintf('WARNING: column %s not found\n', col);
    end
end

% ultraharmonic dB columns
for h = selected_ultraharmonics
    col = sprintf('U%d_db', h);
    if ismember(col, T.Properties.VariableNames)
        col_names{end+1} = col;
        feature_labels{end+1} = sprintf('U%d (%.2f MHz)', h, (h+0.5)*0.5);
    else
        fprintf('WARNING: column %s not found\n', col);
    end
end

% broadband noise
if ismember('broadband_noise_dB', T.Properties.VariableNames)
    col_names{end+1} = 'broadband_noise_dB';
    feature_labels{end+1} = 'BB Noise';
end

n_features = length(col_names);
fprintf('Found %d features:\n', n_features);
for i = 1:n_features
    fprintf('  %d. %s -> %s\n', i, col_names{i}, feature_labels{i});
end

% extract as matrix
X = T{:, col_names};

% remove rows with NaN or Inf
valid = all(isfinite(X), 2);
X = X(valid, :);
n_obs = size(X, 1);
fprintf('Valid samples after cleanup: %d\n', n_obs);

%% create binary label y using 6 dB threshold on broadband noise
bb_col = find(strcmp(feature_labels, 'BB Noise'));
if isempty(bb_col)
    error('BB Noise feature not found — check CSV column names');
end
BB_dB = X(:, bb_col);

baseline_BB = median(BB_dB);
threshold = baseline_BB + 6;

bb_events = double(BB_dB > threshold);
n_events = sum(bb_events);

fprintf('\nBroadband labeling:\n');
fprintf('  Baseline (median): %.2f dB\n', baseline_BB);
fprintf('  Threshold (+6 dB): %.2f dB\n', threshold);
fprintf('  Events: %d / %d (%.1f%%)\n', n_events, n_obs, 100*n_events/n_obs);

% fallback if fewer than 5 events
if n_events < 5
    fprintf('  Too few events — falling back to 75th percentile\n');
    threshold = prctile(BB_dB, 75);
    bb_events = double(BB_dB > threshold);
    n_events = sum(bb_events);
    fprintf('  New threshold: %.2f dB, events: %d\n', threshold, n_events);
end

% shift: features(t) -> label(t+1)
X = X(1:end-1, :);
y = bb_events(2:end);
n_obs = length(y);

fprintf('\nFinal: %d samples x %d features (%d pos, %d neg)\n', ...
    n_obs, n_features, sum(y==1), sum(y==0));

%% train model
mdl = fitrensemble(X, y, 'Method', 'Bag', 'NumLearningCycles', 50);

%% compute SHAP values
fprintf('\nComputing SHAP values for %d samples (this may take a while)...\n', n_obs);
explainer = shapley(mdl, X, 'QueryPoints', X);

% reshape SHAP values
shap_table = explainer.ShapleyValues;
shap_values = shap_table.ShapleyValue;
shap_data = reshape(shap_values, n_obs, n_features);

% normalize SHAP values to [-1, 1]
max_shap = max(abs(shap_data(:)));
shap_data = shap_data / max_shap;

% calculate mean absolute SHAP
mean_abs_shap = mean(abs(shap_data), 1);
[~, sort_idx] = sort(mean_abs_shap, 'descend');

%% plot
set(0, 'DefaultFigureWindowStyle', 'normal');
figure('Position', [100, 100, 1400, 600]);

% left subplot: global feature importance
subplot(1, 2, 1);

hold on;
for i = 1:n_features
    feat_idx = sort_idx(i);
    if feat_idx == bb_col
        bar_color = [0.9, 0.3, 0.3];  % red for BB noise
    elseif startsWith(feature_labels{feat_idx}, 'U')
        bar_color = [0.3, 0.8, 0.4];  % green for ultraharmonics
    else
        bar_color = [0.2, 0.6, 0.8];  % blue for harmonics
    end
    barh(i, mean_abs_shap(feat_idx), 'FaceColor', bar_color);
end
hold off;

set(gca, 'YDir', 'reverse');
title('Global Feature Importance');
xlabel('Mean |SHAP| Value');
ylabel('Feature');
yticks(1:n_features);
yticklabels(feature_labels(sort_idx));
grid on;

% add values as text
for i = 1:n_features
    feat_idx = sort_idx(i);
    x_pos = mean_abs_shap(feat_idx) + 0.005;
    text(x_pos, i, sprintf('%.3f', mean_abs_shap(feat_idx)), ...
        'FontSize', 8, 'VerticalAlignment', 'middle');
end

% right subplot: beeswarm 
subplot(1, 2, 2);
hold on;

for i = 1:n_features
    feat_idx = sort_idx(i);
    feat_shap = shap_data(:, feat_idx);

    jitter = 0.15 * randn(size(feat_shap));

    scatter(feat_shap, i + jitter, 15, 'filled', ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');

    plot([mean(feat_shap), mean(feat_shap)], [i-0.25, i+0.25], ...
        'r-', 'LineWidth', 2);
end

hold off;

set(gca, 'YDir', 'reverse');
title('SHAP Values Distribution');
xlabel('SHAP Value');
ylabel('Feature');
xlim([-1, 1]);
ylim([0.5, n_features + 0.5]);
yticks(1:n_features);
yticklabels(feature_labels(sort_idx));
grid on;
xline(0, 'k--', 'LineWidth', 1);

sgtitle(sprintf('SHAP Analysis: In-vitro Probe Data — %d Features (Reduced)', n_features), 'FontSize', 14);

%% print feature ranking
fprintf('\n=== FEATURE IMPORTANCE RANKING (Dataset 2 — Probe) ===\n');
fprintf('Most Important (Top) → Least Important (Bottom):\n\n');
for i = 1:n_features
    feat_num = sort_idx(i);
    marker = '';
    if feat_num == bb_col
        marker = '  <-- BB Noise';
    end
    fprintf('%d. %s (Mean |SHAP| = %.4f)%s\n', ...
        i, feature_labels{feat_num}, mean_abs_shap(feat_num), marker);
end

%% save fig
print('SHAP_analysis_dataset2_probe.png', '-dpng', '-r300');
