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
    % get first burst data
    first_burst_data = reader.time_series_mV{1};
    
    % check if dimensions are 1 x N (single channel, time points)
    fprintf('first burst data dimensions: %d x %d\n', ...
        size(first_burst_data, 1), size(first_burst_data, 2));
    
    % transpose (275000 x 1 should be 1 x 275000)
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
    
    % create plots (undocked)
    set(0, 'DefaultFigureWindowStyle', 'normal');
    
    % fig 1: single burst time series
    figure('Name', 'time series - first burst', 'Position', [50 50 900 600]);
    
    subplot(2, 2, 1);
    plot(time_ms, first_burst_data);
    title('first burst full signal');
    xlabel('time (ms)');
    ylabel('amplitude (mV)');
    grid on;
    
    subplot(2, 2, 2);
    plot(time_ms(1:min(5000, end)), first_burst_data(1:min(5000, end)));
    title('first 5 ms in detail');
    xlabel('time (ms)');
    ylabel('amplitude (mV)');
    grid on;
    
    subplot(2, 2, 3);
    histogram(first_burst_data, 100);
    title('amplitude distribution');
    xlabel('amplitude (mV)');
    ylabel('count');
    grid on;
    
    subplot(2, 2, 4);
    [pxx, f] = pwelch(first_burst_data, [], [], [], 1/reader.dt);
    plot(f/1e6, 10*log10(pxx));
    title('power spectrum');
    xlabel('frequency (MHz)');
    ylabel('power (dB)');
    grid on;
    xlim([0 10]);
    ylim([-100 0]);
    
    % fig 2: commparing multiple bursts
    figure('Name', 'multiple bursts', 'Position', [1000 50 900 400]);
    
    num_bursts_plot = min(5, reader.n_bursts);
    colors = lines(num_bursts_plot);
    
    hold on;
    for burst = 1:num_bursts_plot
        burst_data = reader.time_series_mV{burst};
        if size(burst_data, 1) > size(burst_data, 2)
            burst_data = burst_data';
        end
        plot(time_ms(1:min(1000, end)), burst_data(1:min(1000, end)), ...
            'Color', colors(burst, :), 'DisplayName', sprintf('burst %d', burst));
    end
    hold off;
    
    title('first 5 bursts (first 1 ms)');
    xlabel('time (ms)');
    ylabel('amplitude (mV)');
    legend('show', 'Location', 'best');
    grid on;
    
    % fig 3: frequency domain data with harmonic analysis
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
    
    % signal analysis
    fprintf('\n=== SIGNAL ANALYSIS ===\n');
    
    % signal statistics
    fprintf('signal statistics:\n');
    fprintf('  mean: %.4f mV\n', mean(first_burst_data));
    fprintf('  standard deviation: %.4f mV\n', std(first_burst_data));
    fprintf('  peak-to-peak: %.4f mV\n', max(first_burst_data) - min(first_burst_data));
    fprintf('  RMS: %.4f mV\n', rms(first_burst_data));
    
    % find ultrasound pulse
    [max_val, max_idx] = max(abs(first_burst_data));
    fprintf('maximum amplitude: %.4f mV at %.3f ms\n', ...
        first_burst_data(max_idx), time_ms(max_idx));
    
    % calculate SNR (crude estimate)
    noise_std = std(first_burst_data(1:round(0.1*length(first_burst_data))));
    signal_std = std(first_burst_data);
    snr_estimate = 20*log10(signal_std/noise_std);
    fprintf('estimated SNR: %.2f dB\n', snr_estimate);
    
else
    fprintf('oop no time series data found.\n');
end