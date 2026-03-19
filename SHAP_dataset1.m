%% PCD SHAP Analysis
% training a model to predict broadband emission from harmonic features
% SHAP shows which features push toward/away from broadband (inertial) events
% positive SHAP = pushes toward inertial (unsafe)
% negative SHAP = pushes toward safe
clear; clc; close all;
rng(42);

load('PCD_Features_AllFiles.mat', 'Results');

% build arrays
H_dB  = vertcat(Results.harmonic_power_dB);        % N x 4
U_dB  = vertcat(Results.ultraharmonic_power_dB);   % N x 3
BB_dB = [Results.broadband_noise_power_dB]';       % N x 1
N = length(BB_dB);

% detect sweeps and normalize BB for labeling
% detect concentration sweep boundaries by looking for large drops in 
% f0 power, and normalize BB within each sweep relative to that sweep's
% baseline (diff concentrations might have diff baseline noise floors)
f0_dB = H_dB(:, 1);
f0_diff = diff(f0_dB);
reset_indices = find(f0_diff < -10);
sweep_starts = [1; reset_indices + 1];
sweep_ends   = [reset_indices; N];
n_sweeps = length(sweep_starts);
n_baseline_per_sweep = 3;

BB_normalized = zeros(N, 1);
for s = 1:n_sweeps
    idx = sweep_starts(s):sweep_ends(s);
    base_idx = sweep_starts(s):min(sweep_starts(s) + n_baseline_per_sweep - 1, sweep_ends(s));
    sweep_baseline = mean(BB_dB(base_idx));
    BB_normalized(idx) = BB_dB(idx) - sweep_baseline;
end

%% MODEL SETUP / feature matrix preparation
%  input features: H1-H3, U1-U3 (NO broadband — we're predicting it)
%  target: broadband noise elevation 

X = [H_dB(:, 2:4), U_dB];

feature_names = {'H1 (1.0 MHz)', 'H2 (1.5 MHz)', 'H3 (2.0 MHz)', ...
                 'U1 (0.75 MHz)', 'U2 (1.25 MHz)', 'U3 (1.75 MHz)'};
n_features = length(feature_names);

% remove NaN
valid = all(~isnan(X), 2);
X = X(valid, :);
y = BB_normalized(valid);  % continuous target: dB above baseline
n_obs = sum(valid);

fprintf('Samples: %d, Features: %d\n', n_obs, n_features);
fprintf('Target range: %.2f to %.2f dB above baseline\n', min(y), max(y));

%% train regression ensemble (predicting continuous BB elevation)
mdl = fitrensemble(X, y, 'Method', 'Bag', 'NumLearningCycles', 100);

% evaluate fit
y_pred = predict(mdl, X);
r2 = 1 - sum((y - y_pred).^2) / sum((y - mean(y)).^2);
rmse = sqrt(mean((y - y_pred).^2));
fprintf('R² = %.3f, RMSE = %.3f dB\n', r2, rmse);

%% compute SHAP values
fprintf('Computing SHAP values...\n');
explainer = shapley(mdl, X, 'QueryPoints', X);
shap_table = explainer.ShapleyValues;
shap_values = reshape(shap_table.ShapleyValue, n_obs, n_features);

% mean absolute SHAP for global importance ranking (on normalized values)
mean_abs_shap = mean(abs(shap_values), 1);
[~, sort_idx] = sort(mean_abs_shap, 'descend');

% mean signed SHAP for directionality
mean_signed_shap = mean(shap_values, 1);

%% plot global feature importance + beeswarm plot
figure('Position', [100 100 1400 550]);

% left plot = global importance (absolute SHAP) 
subplot(1, 2, 1);
hold on;
for i = 1:n_features
    fi = sort_idx(i);
    barh(i, mean_abs_shap(fi), 'FaceColor', [0.2 0.6 0.8]);
end
set(gca, 'YDir', 'reverse');
title('Global Feature Importance', 'FontSize', 16);
xlabel('Mean |SHAP| Value', 'FontSize', 14);
ylabel('Feature', 'FontSize', 14);
yticks(1:n_features); yticklabels(feature_names(sort_idx));
grid on; set(gca, 'FontSize', 14);

for i = 1:n_features
    fi = sort_idx(i);
    text(mean_abs_shap(fi) + 0.005 * max(mean_abs_shap), i, ...
        sprintf('%.4f', mean_abs_shap(fi)), ...
        'FontSize', 12, 'VerticalAlignment', 'middle');
end

% right plot = beeswarm colored by feature value 
% red = high feature value, blue = low
subplot(1, 2, 2);
hold on;
for i = 1:n_features
    feat_idx = sort_idx(i);
    feat_shap = shap_values(:, feat_idx);
    feat_vals = X(:, feat_idx);

    % normalize feature values to [0, 1] for coloring
    feat_norm = (feat_vals - min(feat_vals)) / (max(feat_vals) - min(feat_vals) + eps);

    jitter = 0.15 * randn(n_obs, 1);
    scatter(feat_shap, i + jitter, 20, feat_norm, 'filled', 'MarkerFaceAlpha', 0.6);
end

colormap(turbo);
cb = colorbar;
cb.Label.String = 'Feature Value (Low → High)';
cb.Label.FontSize = 13;

set(gca, 'YDir', 'reverse');
xline(0, 'k--', 'LineWidth', 1.5);
title('SHAP Values (→ = more broadband/inertial)', 'FontSize', 16);
xlabel('SHAP Value (← Safe | Unsafe →)', 'FontSize', 14);
yticks(1:n_features); yticklabels(feature_names(sort_idx));
xlim([-1, 1]);
ylim([0.5, n_features + 0.5]);
grid on; set(gca, 'FontSize', 14);

sgtitle(sprintf('PCD SHAP Analysis — Broadband Prediction (R² = %.3f)', r2), 'FontSize', 16);
print('SHAP_study_replication_summary.png', '-dpng', '-r300');
%% print results
fprintf('\n=== FEATURE RANKING ===\n');
fprintf('Positive SHAP = pushes toward higher broadband (inertial/unsafe)\n');
fprintf('Negative SHAP = pushes toward lower broadband (safe)\n\n');
fprintf('%-5s | %-20s | |SHAP|  | Signed   | Interpretation\n', 'Rank', 'Feature');
fprintf('------|----------------------|---------|----------|------------------\n');
for i = 1:n_features
    fi = sort_idx(i);
    if mean_signed_shap(fi) > 0.001
        interp = 'Higher value → MORE broadband';
    elseif mean_signed_shap(fi) < -0.001
        interp = 'Higher value → LESS broadband';
    else
        interp = 'Mixed / nonlinear relationship';
    end
    fprintf('  %d   | %-20s | %.4f  | %+.4f  | %s\n', ...
        i, feature_names{fi}, mean_abs_shap(fi), mean_signed_shap(fi), interp);
end