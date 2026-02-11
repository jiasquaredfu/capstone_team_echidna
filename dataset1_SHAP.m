clear; clc; close all;

% load features extracted from dataset1_pcd_allfiles.m
load('PCD_features_for_SHAP.mat');
X_real = X(:, 1:6);

% train model
mdl = fitrensemble(X_real, y, 'Method', 'Bag', 'NumLearningCycles', 50);

% compute SHAP values
explainer = shapley(mdl, X_real, 'QueryPoints', X_real);

% get SHAP values table
shap_table = explainer.ShapleyValues;

% table is in long format with columns: 'Predictor' and 'ShapleyValue'
% reshape to wide format
predictors = shap_table.Predictor;
shap_values = shap_table.ShapleyValue;

% reshape to matrix: each column is a predictor, each row is an observation
n_obs = length(y);
shap_data = reshape(shap_values, n_obs, 6);

% calculate mean absolute SHAP for each feature
mean_abs_shap = mean(abs(shap_data), 1);

% sort features by importance (descending: most important first)
[~, sort_idx] = sort(mean_abs_shap, 'descend');

% define feature names with frequencies
feature_names = {'H1 (1.0 MHz)', 'H2 (1.5 MHz)', 'H3 (2.0 MHz)', ...
                 'U1.5 (0.75 MHz)', 'U2.5 (1.25 MHz)', 'U3.5 (1.75 MHz)'};

% plot
figure('Position', [100, 100, 1200, 500]);

% left subplot: global feature importance
subplot(1, 2, 1);
barh(mean_abs_shap(sort_idx), 'FaceColor', [0.2, 0.6, 0.8]);
set(gca, 'YDir', 'reverse');  % put most important feature at the top
title('Global Feature Importance');
xlabel('Mean |SHAP| Value');
ylabel('Feature');
yticks(1:6);
% Use the actual feature names sorted by importance
yticklabels(feature_names(sort_idx));
grid on;

% right subplot: swarm plot of SHAP values
subplot(1, 2, 2);
hold on;
for i = 1:6
    feat_idx = sort_idx(i);  % current feature in ranking
    feat_shap = shap_data(:, feat_idx);
    
    % add jitter
    jitter = 0.15 * randn(size(feat_shap));
    
    % plot points
    scatter(feat_shap, i + jitter, 25, 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k');
    
    % add mean line
    plot([mean(feat_shap), mean(feat_shap)], [i-0.25, i+0.25], ...
        'r-', 'LineWidth', 2);
end
hold off;

set(gca, 'YDir', 'reverse');  % set most important feature at top
title('SHAP Values Distribution');
xlabel('SHAP Value');
ylabel('Feature');
ylim([0.5, 6.5]);
yticks(1:6);

% use the same feature names sorted by importance
yticklabels(feature_names(sort_idx));
grid on;
xline(0, 'k--', 'LineWidth', 1);  % Add zero line

sgtitle('SHAP Analysis: In-vitro PCD');

% add mean SHAP values as text on the bar plot
ax1 = subplot(1, 2, 1);
for i = 1:6
    feat_idx = sort_idx(i);
    x_pos = mean_abs_shap(feat_idx) + 0.1;  % Slight offset from bar
    text(x_pos, i, sprintf('%.2f', mean_abs_shap(feat_idx)), ...
        'FontSize', 9, 'VerticalAlignment', 'middle');
end

fprintf('\n=== FEATURE IMPORTANCE RANKING ===\n');
fprintf('Most Important (Top) → Least Important (Bottom):\n\n');
for i = 1:6
    rank = i;
    feat_num = sort_idx(i);
    fprintf('%d. %s (Mean |SHAP| = %.4f)\n', ...
        rank, feature_names{feat_num}, mean_abs_shap(feat_num));
end

% save figure with high resolution
print('SHAP_analysis_dataset1.png', '-dpng', '-r300');