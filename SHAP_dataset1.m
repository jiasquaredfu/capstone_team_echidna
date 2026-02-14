clear; clc; close all;
rng(42);  % fix random seed for reproducibility

% load features directly from extraction script output
load('PCD_Features_AllFiles.mat', 'Results');

% Build arrays from Results struct
H_dB  = vertcat(Results.harmonic_power_dB);        % N x 4 (f0, H1, H2, H3)
U_dB  = vertcat(Results.ultraharmonic_power_dB);   % N x 3
BB_dB = [Results.broadband_noise_power_dB]';       % N x 1

% feature matrix: drop f0 (col 1), keep H1-H3 (cols 2:4), add U's and BB
%   Col 1: H1 (1.0 MHz)
%   Col 2: H2 (1.5 MHz)
%   Col 3: H3 (2.0 MHz)
%   Col 4: U1 (0.75 MHz)
%   Col 5: U2 (1.25 MHz)
%   Col 6: U3 (1.75 MHz)
%   Col 7: Broadband Noise  
X = [H_dB(:, 2:4), U_dB, BB_dB];

% remove rows with NaN
valid = all(~isnan(X), 2);
X = X(valid, :);
BB_valid = BB_dB(valid);

%% create binary label y using paper's 6 dB threshold
%  shifted: features at t predict broadband event at t+1
baseline_BB = median(BB_valid);
threshold = baseline_BB + 6;

bb_events = double(BB_valid > threshold);

fprintf('Broadband labeling: baseline=%.2f dB, threshold=%.2f dB\n', ...
    baseline_BB, threshold);
fprintf('  Events: %d / %d (%.1f%%)\n', sum(bb_events), length(bb_events), ...
    100*sum(bb_events)/length(bb_events));

% fallback if fewer than 5 events
if sum(bb_events) < 5
    fprintf('  Too few events — falling back to 75th percentile\n');
    threshold = prctile(BB_valid, 75);
    bb_events = double(BB_valid > threshold);
    fprintf('  New threshold: %.2f dB, events: %d\n', threshold, sum(bb_events));
end

% shift: features(t) -> label(t+1)
X = X(1:end-1, :);
y = bb_events(2:end);

n_features = size(X, 2);  % 7
n_obs = length(y);

fprintf('Final: %d samples x %d features (%d pos, %d neg)\n', ...
    n_obs, n_features, sum(y==1), sum(y==0));

%% train model
mdl = fitrensemble(X, y, 'Method', 'Bag', 'NumLearningCycles', 50);

%% compute SHAP values
explainer = shapley(mdl, X, 'QueryPoints', X);

% get SHAP values and reshape
shap_table = explainer.ShapleyValues;
shap_values = shap_table.ShapleyValue;
shap_data = reshape(shap_values, n_obs, n_features);

% normalize SHAP values to [-1, 1]
max_shap = max(abs(shap_data(:)));
shap_data = shap_data / max_shap;

% calculate mean absolute SHAP for each feature
mean_abs_shap = mean(abs(shap_data), 1);
[~, sort_idx] = sort(mean_abs_shap, 'descend');

%% define feature names
feature_names = {'H1 (1.0 MHz)', 'H2 (1.5 MHz)', 'H3 (2.0 MHz)', ...
                 'U1 (0.75 MHz)', 'U2 (1.25 MHz)', 'U3 (1.75 MHz)', ...
                 'BB Noise'};

%% plot
set(0, 'DefaultFigureWindowStyle', 'normal');  % undock
figure('Position', [100, 100, 1300, 550]);

% left subplot: global feature importance 
subplot(1, 2, 1);

% red for BB noise, blue for everything else
bar_colors = repmat([0.2, 0.6, 0.8], n_features, 1);
for i = 1:n_features
    if sort_idx(i) == 7  % BB noise index
        bar_colors(i, :) = [0.9, 0.3, 0.3];
    end
end

% plot bar colors individually
hold on;
for i = 1:n_features
    barh(i, mean_abs_shap(sort_idx(i)), 'FaceColor', bar_colors(i,:));
end
hold off;

set(gca, 'YDir', 'reverse');
title('Global Feature Importance');
xlabel('Mean |SHAP| Value');
ylabel('Feature');
yticks(1:n_features);
yticklabels(feature_names(sort_idx));
grid on;

% right subplot: swarm plot 
subplot(1, 2, 2);
hold on;
for i = 1:n_features
    feat_idx = sort_idx(i);
    feat_shap = shap_data(:, feat_idx);

    % jitter
    jitter = 0.15 * randn(size(feat_shap));

    % plot points
    scatter(feat_shap, i + jitter, 25, 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');

    % add mean line
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
yticklabels(feature_names(sort_idx));
grid on;
xline(0, 'k--', 'LineWidth', 1);

sgtitle('SHAP Analysis: In-vitro PCD (with Broadband Noise)');

% add mean SHAP values as text on bar plot
subplot(1, 2, 1);
for i = 1:n_features
    feat_idx = sort_idx(i);
    x_pos = mean_abs_shap(feat_idx) + 0.02 * max(mean_abs_shap);
    text(x_pos, i, sprintf('%.2f', mean_abs_shap(feat_idx)), ...
        'FontSize', 9, 'VerticalAlignment', 'middle');
end

%% print feature ranking
fprintf('\n=== FEATURE IMPORTANCE RANKING ===\n');
fprintf('Most Important (Top) → Least Important (Bottom):\n\n');
for i = 1:n_features
    feat_num = sort_idx(i);
    marker = '';
    if feat_num == 7
        marker = '  <-- BB Noise';
    end
    fprintf('%d. %s (Mean |SHAP| = %.4f)%s\n', ...
        i, feature_names{feat_num}, mean_abs_shap(feat_num), marker);
end

%% save fig
print('SHAP_analysis_dataset1.png', '-dpng', '-r300');