%% PCD Time-Frequency Analysis — H, U, BB vs Acoustic Pressure
% plots total harmonic, ultraharmonic, and broadband noise power
% as 3 separate lines across the pressure ramp
% one figure per concentration folder
clear; clc; close all;
set(0, 'DefaultFigureWindowStyle', 'normal');

pcd_base     = '/Users/stephanieliu/Desktop/in vitro PCD & probe/pcd';
ambient_path = fullfile(pcd_base, 'ambient');

data_folders = { ...
    fullfile(pcd_base, '100x_01');   ...
    fullfile(pcd_base, '100x_02');   ...
    fullfile(pcd_base, '1000x_01');  ...
    fullfile(pcd_base, '1000x_02') };

% parameters
drive_freq   = 0.5e6;
n_harmonics  = 4;
uh_orders    = [1.5 2.5 3.5];
bin_range    = 1000;
noise_exclusion = 50;
uh_upper_harm   = [2, 3, 4];  % BB band: U_k to H_(uh_upper_harm(k))

pressures = 0.025 : 0.025 : 0.800;  % 32 levels (MPa)
n_pres = length(pressures);

% ambient files
amb_files = dir(fullfile(ambient_path, '*.mat'));
amb_files_alt = amb_files(2:2:end);

% process each folder
for d = 1:length(data_folders)
    [~, folder_label] = fileparts(data_folders{d});
    fprintf('\nProcessing %s...\n', folder_label);

    mat_files = dir(fullfile(data_folders{d}, '*.mat'));
    mat_files = mat_files(2:2:end);
    n_files = min(length(mat_files), n_pres);

    % storage: one value per pressure
    H_total_dB  = NaN(n_files, 1);  % sum of all harmonic power
    U_total_dB  = NaN(n_files, 1);  % sum of all ultraharmonic power
    BB_total_dB = NaN(n_files, 1);  % average broadband noise power

    for i = 1:n_files
        reader = load(fullfile(data_folders{d}, mat_files(i).name));
        fs = 1 / reader.Tinterval;
        [pxx, f] = pwelch(reader.B, [], [], [], fs);
        df = f(2) - f(1);
        spec = sqrt(pxx);

        % harmonics: total power across f0, H1, H2, H3 
        h_power_total = 0;
        harm_idx = zeros(1, n_harmonics);
        for k = 1:n_harmonics
            [~, harm_idx(k)] = min(abs(f - k * drive_freq));
            startBin = max(1, harm_idx(k) - bin_range);
            endBin   = min(length(pxx), harm_idx(k) + bin_range);
            h_power_total = h_power_total + sum(pxx(startBin:endBin)) * df;
        end
        H_total_dB(i) = 10*log10(h_power_total + eps);

        % ultraharmonics: total power across U1.5, U2.5, U3.5 
        u_power_total = 0;
        search_range = bin_range * df;
        uh_idx = zeros(1, length(uh_orders));
        for k = 1:length(uh_orders)
            f_target = uh_orders(k) * drive_freq;
            idx_win = find(f >= f_target - search_range & f <= f_target + search_range);
            if isempty(idx_win), continue; end
            [~, max_idx] = max(spec(idx_win));
            uh_idx(k) = idx_win(max_idx);
            startBin = max(1, uh_idx(k) - bin_range);
            endBin   = min(length(pxx), uh_idx(k) + bin_range);
            u_power_total = u_power_total + sum(pxx(startBin:endBin)) * df;
        end
        U_total_dB(i) = 10*log10(u_power_total + eps);

        % broadband noise: average across all 3 BB bands 
        bb_powers = [];
        for k = 1:length(uh_orders)
            if uh_idx(k) == 0, continue; end
            bb_start = max(1, uh_idx(k)+ noise_exclusion + 1);
            bb_end   = min(length(pxx), harm_idx(uh_upper_harm(k)) - noise_exclusion - 1);
            if bb_end <= bb_start, continue; end
            bb_powers(end+1) = mean(pxx(bb_start:bb_end)) * df;
        end
        if ~isempty(bb_powers)
            BB_total_dB(i) = 10*log10(mean(bb_powers) + eps);
        end
    end

    % plot
    p = pressures(1:n_files);

    figure('Color', 'w', 'Position', [100 + d*50, 100 + d*50, 800, 500]);
    hold on;

    plot(p, H_total_dB,  'r-o', 'LineWidth', 2, 'MarkerSize', 5, ...
        'MarkerFaceColor', 'r', 'DisplayName', 'Harmonic (H)');
    plot(p, U_total_dB,  'g-s', 'LineWidth', 2, 'MarkerSize', 5, ...
        'MarkerFaceColor', [0 0.5 0], 'DisplayName', 'Ultraharmonic (U)');
    plot(p, BB_total_dB, 'b-^', 'LineWidth', 2, 'MarkerSize', 5, ...
        'MarkerFaceColor', 'b', 'DisplayName', 'Broadband Noise (BB)');

    xlabel('Acoustic Pressure (MPa)', 'FontSize', 25, 'FontWeight', 'bold');
    ylabel('Total Power (dB)', 'FontSize', 24, 'FontWeight', 'bold');
    title(sprintf('Time-Frequency Analysis — %s', folder_label), ...
        'FontSize', 24, 'FontWeight', 'bold', 'Interpreter', 'none');
    legend('Location', 'northwest', 'FontSize', 11);
    grid on;
    set(gca, 'FontSize', 24, 'LineWidth', 1.2);
    hold off;
end
