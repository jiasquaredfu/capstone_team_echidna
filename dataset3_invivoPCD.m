%% in vivo PCD dataset analysis
clear; clc; close all;

fprintf('current folder: %s\n', pwd);

% add repo to pat
repo_path = '04_Projects/17_FUS_Instruments_Interface';
if exist(repo_path, 'dir')
    addpath(genpath(repo_path));
    fprintf('added repository to path\n');
else
    error('repo not found at: %s', repo_path);
end

% check RKPC_data_reader
if exist('RKPC_data_reader.m', 'file')
    fprintf('RKPC_data_reader.m found\n');
else
    error('RKPC_data_reader.m not found.');
end

% find HDF5 files in current folder
h5_files = dir('*.hdf5');
if isempty(h5_files)
    h5_files = dir('*.h5');
end

if isempty(h5_files)
    error('no HDF5 files found in current folder.');
end

fprintf('\nfound %d HDF5 files:\n', length(h5_files));
for i = 1:length(h5_files)
    fprintf('%d. %s\n', i, h5_files(i).name);
end
%% look at the first .hdf5 file
fprintf('\nloading file: %s\n', h5_files(1).name);

% call reader
reader = RKPC_data_reader(h5_files(1).name);

% data is now stored in the reader object's properties
fprintf('data loaded successfully into reader object\n');

% display basic information from the reader
fprintf('\n=== DATA INFORMATION ===\n');
fprintf('file name: %s\n', reader.file_name);
fprintf('number of bursts: %d\n', reader.n_bursts);
fprintf('number of bursts planned: %d\n', reader.n_bursts_planned);
fprintf('drive frequency: %.2f MHz\n', reader.drive_frequency_MHz);
fprintf('start pressure: %.2f MPa\n', reader.start_pressure_MPa);
fprintf('burst length: %.2f ms\n', reader.burst_length_ms);
fprintf('burst period: %.2f ms\n', reader.burst_period_ms);
fprintf('sampling frequency: %.2f MHz\n', reader.fs_MHz);
fprintf('PRF: %.2f Hz\n', reader.prf_Hz);
fprintf('duty cycle: %.4f\n', reader.duty_cycle);
fprintf('time step (dt): %.2e seconds\n', reader.dt);
fprintf('frequency step (df): %.2f MHz\n', reader.df_MHz);

% check data dimensions
fprintf('\n=== TIME SERIES DATA ===\n');
if reader.n_bursts > 0 && ~isempty(reader.time_series_mV{1})
    % get first burst data for reference dimensions
    first_burst_data = reader.time_series_mV{1};
    
    % check if dimensions are 1 x N (single channel, time points)
    fprintf('first burst data dimensions: %d x %d\n', ...
        size(first_burst_data, 1), size(first_burst_data, 2));
    
    % transpose if needed
    if size(first_burst_data, 1) > size(first_burst_data, 2)
        fprintf('transposing data (single channel, time points)\n');
        first_burst_data = first_burst_data';
        fprintf('after transpose: %d x %d\n', ...
            size(first_burst_data, 1), size(first_burst_data, 2));
    end
    
    fprintf('number of time points: %d\n', size(first_burst_data, 2));
    fprintf('time duration: %.3f ms\n', size(first_burst_data, 2) * reader.dt * 1000);
    
    % create time vector in milliseconds
    time_ms = (0:size(first_burst_data, 2)-1) * reader.dt * 1000;
    
    % collect all bursts into matrix (aggregate)
    n_bursts = reader.n_bursts;
    all_bursts_matrix = zeros(n_bursts, length(first_burst_data));
    
    for burst = 1:n_bursts
        burst_data = reader.time_series_mV{burst};
        if size(burst_data, 1) > size(burst_data, 2)
            burst_data = burst_data';
        end
        all_bursts_matrix(burst, :) = burst_data;
    end
    
    % calculate aggregated statistics
    mean_signal = mean(all_bursts_matrix, 1);
    std_signal = std(all_bursts_matrix, 0, 1);
    median_signal = median(all_bursts_matrix, 1);
    
    fprintf('\naggregation statistics:\n');
    fprintf('  averaging across %d bursts\n', n_bursts);
    fprintf('  single channel (no spatial channels)\n');
    
    % create plots (undocked)
    set(0, 'DefaultFigureWindowStyle', 'normal');
    
    % fig 1: aggregated time series analysis
    figure('Name', 'aggregated time series analysis', 'Position', [50 50 1200 800]);
    
    % plot 1: mean signal ± std
    subplot(2, 3, 1);
    plot(time_ms, mean_signal, 'b-', 'LineWidth', 1.5);
    hold on;
    plot(time_ms, mean_signal + std_signal, 'r--', 'LineWidth', 0.5);
    plot(time_ms, mean_signal - std_signal, 'r--', 'LineWidth', 0.5);
    fill([time_ms, fliplr(time_ms)], ...
         [mean_signal + std_signal, fliplr(mean_signal - std_signal)], ...
         'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold off;
    title(sprintf('mean signal ± std (n=%d bursts)', n_bursts));
    xlabel('time (ms)');
    ylabel('amplitude (mV)');
    grid on;
    legend('mean', '±1 std', 'Location', 'best');
    
    % plot 2: first 5 ms detail
    subplot(2, 3, 2);
    plot(time_ms(1:min(5000, end)), mean_signal(1:min(5000, end)), 'b-', 'LineWidth', 1.5);
    hold on;
    fill([time_ms(1:min(5000, end)), fliplr(time_ms(1:min(5000, end)))], ...
         [mean_signal(1:min(5000, end)) + std_signal(1:min(5000, end)), ...
          fliplr(mean_signal(1:min(5000, end)) - std_signal(1:min(5000, end)))], ...
         'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold off;
    title('first 5 ms detail');
    xlabel('time (ms)');
    ylabel('amplitude (mV)');
    grid on;
    
    % plot 3: burst-to-burst variability heatmap
    subplot(2, 3, 3);
    imagesc(time_ms, 1:n_bursts, all_bursts_matrix);
    title('burst heatmap (burst # vs time)');
    xlabel('time (ms)');
    ylabel('burst number');
    colorbar;
    ylabel(colorbar, 'amplitude (mV)');
    
    % plot 4: amplitude distribution of aggregated signal
    subplot(2, 3, 4);
    histogram(mean_signal, 50);
    title('amplitude distribution (mean signal)');
    xlabel('amplitude (mV)');
    ylabel('count');
    grid on;
    
    % plot 5: power spectrum of mean signal
    subplot(2, 3, 5);
    [pxx_mean, f] = pwelch(mean_signal, [], [], [], 1/reader.dt);
    plot(f/1e6, 10*log10(pxx_mean), 'b-', 'LineWidth', 1.5);
    title('power spectrum (mean signal)');
    xlabel('frequency (MHz)');
    ylabel('power (dB)');
    grid on;
    xlim([0 10]);
    ylim([-100 0]);
    
    % plot 6: burst statistics over time
    subplot(2, 3, 6);
    burst_peaks = max(abs(all_bursts_matrix), [], 2);
    burst_rms = rms(all_bursts_matrix, 2);
    
    plot((1:n_bursts)*reader.burst_period_ms/1000, burst_peaks, 'r-o', 'MarkerSize', 4);
    hold on;
    plot((1:n_bursts)*reader.burst_period_ms/1000, burst_rms, 'b-s', 'MarkerSize', 4);
    hold off;
    title('burst amplitude over experiment');
    xlabel('experiment time (s)');
    ylabel('amplitude (mV)');
    legend('peak amplitude', 'RMS amplitude', 'Location', 'best');
    grid on;
    
    % fig 2: burst comparsion
    figure('Name', 'burst-to-burst comparison', 'Position', [1000 50 1200 600]);
    
    % select representative bursts
    burst_indices = [1, round(n_bursts/4), round(n_bursts/2), round(3*n_bursts/4), n_bursts];
    burst_indices = burst_indices(burst_indices <= n_bursts);
    colors = lines(length(burst_indices));
    
    % plot 1: selected bursts overlaid
    subplot(1, 2, 1);
    hold on;
    for i = 1:length(burst_indices)
        burst_idx = burst_indices(i);
        plot(time_ms(1:min(5000, end)), ...
             all_bursts_matrix(burst_idx, 1:min(5000, end)), ...
             'Color', colors(i, :), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('burst %d', burst_idx));
    end
    hold off;
    title('selected bursts (first 5 ms)');
    xlabel('time (ms)');
    ylabel('amplitude (mV)');
    legend('show', 'Location', 'best');
    grid on;
    
    % plot 2: correlation between early and late bursts
    subplot(1, 2, 2);
    early_burst = all_bursts_matrix(1, :);
    late_burst = all_bursts_matrix(end, :);
    
    plot(early_burst(1:min(5000, end)), late_burst(1:min(5000, end)), '.');
    hold on;
    plot([-100 100], [-100 100], 'r--', 'LineWidth', 1); % identity line
    hold off;
    title(sprintf('burst correlation: 1 vs %d', n_bursts));
    xlabel('burst 1 amplitude (mV)');
    ylabel(sprintf('burst %d amplitude (mV)', n_bursts));
    grid on;
    axis equal;
    xlim([-100 100]);
    ylim([-100 100]);
    
    % calculate correlation coefficient
    corr_coeff = corrcoef(early_burst, late_burst);
    text(-90, 80, sprintf('r = %.3f', corr_coeff(1,2)), ...
         'FontSize', 12, 'FontWeight', 'bold');
    
    % fig 3: frequency domain with features labeled
    figure('Name', 'cavitation spectrum analysis', 'Position', [50 700 1200 500]);
    
    if ~isempty(reader.freq_spec_dB{1})
        freq_data = reader.freq_spec_dB{1};
        freq_axis = reader.freq_axis;
        f0 = reader.drive_frequency_MHz; % 0.5 MHz
        
        % check dimensions
        if size(freq_data, 1) > size(freq_data, 2)
            freq_data = freq_data';
        end
        
        % plot 1: linear frequency scale with markers
        subplot(1, 2, 1);
        plot(freq_axis, freq_data, 'b-', 'LineWidth', 1.5);
        hold on;
        
        % mark harmonics (H1, H2, H3, H4, H5)
        harmonics = 1:5;
        for h = harmonics
            h_freq = h * f0;
            [~, idx] = min(abs(freq_axis - h_freq));
            plot(h_freq, freq_data(idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
            text(h_freq, freq_data(idx)+5, sprintf('H%d', h), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
        
        % mark ultraharmonics (U1, U2, U3, U4, U5)
        % ultraharmonics are at (n + 1/2) * f0
        ultraharmonics = 1:5;
        for u = ultraharmonics
            u_freq = (u + 0.5) * f0;
            [~, idx] = min(abs(freq_axis - u_freq));
            plot(u_freq, freq_data(idx), 'gs', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
            text(u_freq, freq_data(idx)+5, sprintf('U%d', u), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
        
        % shade broadband noise region 
        % (typically 3-10 MHz for microbubbles)
        broadband_start = 3; % MHz
        broadband_end = 10; % MHz
        y_limits = ylim;
        fill([broadband_start broadband_end broadband_end broadband_start], ...
             [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
             [0.9 0.9 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
        hold off;
        title('cavitation spectrum with harmonics');
        xlabel('frequency (MHz)');
        ylabel('magnitude (dB)');
        grid on;
        xlim([0 10]);
        legend('spectrum', 'harmonics', 'ultraharmonics', 'Location', 'best');
        
        % plot 2: normalized frequency (divide by f0)
        subplot(1, 2, 2);
        plot(freq_axis/f0, freq_data, 'b-', 'LineWidth', 1.5);
        hold on;
        
        % mark harmonics on normalized scale
        for h = harmonics
            plot(h, interp1(freq_axis/f0, freq_data, h), 'ro', ...
                'MarkerSize', 10, 'MarkerFaceColor', 'r');
            text(h, interp1(freq_axis/f0, freq_data, h)+5, sprintf('H%d', h), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
        
        for u = ultraharmonics
            u_norm = u + 0.5;
            plot(u_norm, interp1(freq_axis/f0, freq_data, u_norm), 'gs', ...
                'MarkerSize', 10, 'MarkerFaceColor', 'g');
            text(u_norm, interp1(freq_axis/f0, freq_data, u_norm)+5, sprintf('U%d', u), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
        
        hold off;
        title('normalized frequency (f/f_0)');
        xlabel('frequency / f_0 (f_0 = 0.5 MHz)');
        ylabel('magnitude (dB)');
        grid on;
        xlim([0 20]);
        
    end
    
    % signal analysis using aggregated burst data
    fprintf('\n=== AGGREGATED SIGNAL ANALYSIS ===\n');
    
    % signal statistics on mean signal
    fprintf('mean signal statistics (aggregated across %d bursts):\n', n_bursts);
    fprintf('  mean: %.4f mV\n', mean(mean_signal));
    fprintf('  standard deviation: %.4f mV\n', std(mean_signal));
    fprintf('  peak-to-peak: %.4f mV\n', max(mean_signal) - min(mean_signal));
    fprintf('  RMS: %.4f mV\n', rms(mean_signal));
    
    % find ultrasound pulse in mean signal
    [max_val, max_idx] = max(abs(mean_signal));
    fprintf('maximum amplitude: %.4f mV at %.3f ms\n', ...
        mean_signal(max_idx), time_ms(max_idx));
    
    % calculate SNR using aggregated signal
    noise_std = std(mean_signal(1:round(0.1*length(mean_signal))));
    signal_std = std(mean_signal);
    snr_estimate = 20*log10(signal_std/noise_std);
    fprintf('estimated SNR (aggregated): %.2f dB\n', snr_estimate);
    
    % burst-to-burst variability statistics
    fprintf('\nburst-to-burst variability:\n');
    fprintf('  peak amplitude mean ± std: %.2f ± %.2f mV\n', ...
        mean(burst_peaks), std(burst_peaks));
    fprintf('  RMS amplitude mean ± std: %.2f ± %.2f mV\n', ...
        mean(burst_rms), std(burst_rms));
    fprintf('  coefficient of variation (peak): %.2f%%\n', ...
        100*std(burst_peaks)/mean(burst_peaks));
    
else
    fprintf('no time series data found.\n');
end