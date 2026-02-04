%% Dataset 2 Analysis 
% Capstone II - Team E(Chidna)
% JiaJia Fu 
% 1/25/25

clear;
clc;
close all;


% Add path of imaging function 
function_path = '/Users/jiasquared/Desktop/CODING/capstone_team_echidna/04_Projects/';
addpath(function_path);

% Add path of Philips-ATL L12-5 Probe Data to read 
probe_path = '/Users/jiasquared/Desktop/CODING/capstone_team_echidna/1_and_2_2026-01-16_single_tube_sweep/probe';


%% Extract processsed data for 1st datafile

%all_data(k).rf_data(:, :, frame_idx) = double(data.data{frame_idx});


function out = extract_probe_data(probe_path)
 % Extracts and converts aquisition_reader data to double arrays
        %
        % Input: 
        %       - probe_path: path to probe data for
        %       imaging.common.acquisition_reader to read 
        % 
        % Output:
        %       - out.rf_data: A cell array with each frame of data in
        %       double array form
        %       - out.t [N x 1] time vector 
        %       - out.frame_t [N x 1] timestamps
        %       - out.fs sampling frequency (Hz)

        % Get names of .pacq files in folder 
        data_files = dir(fullfile(probe_path, '*.pacq*'));

        % Apply aquisition_reader function to get data 
        file_path = string(fullfile(data_files(1).folder, data_files(1).name));
        acq = imaging.common.acquisition_reader(file_path);
        disp(acq.data_params)
        disp(acq)
        
        % Define dimensions for array
        N_t = acq.rfdata_dimensions(1);
        N_ch = acq.rfdata_dimensions(2);
        N_frames = acq.number_of_frames;
    
        % Empty array to store rf_data
        rf_data = zeros(N_t, N_ch, N_frames, 'double');
    
        % Convert each frame
        for k = 1:N_frames
            rf_data(:, :, k) = double(acq.data{k});
        end
    
        % Time axis for samples
        fs = acq.data_params.fs;
        t = (0:N_t-1).' / fs;
    
        % Frame timestamps (relative, seconds)
        frame_t = acq.elapsed_timevector(:);
    
        % Pack output
        out = struct();
        out.rf_data = rf_data;
        out.t = t;
        out.frame_t = frame_t;
        out.fs = fs;
end

clean_data = extract_probe_data(probe_path)
fs = clean_data.fs


%% plot raw fs vs time 

frame_1 = clean_data.rf_data(:, 64, 1); % middle channels frame 1 sample 
%{
n_samples = length(frame_1);
t = (0:n_samples-1) / fs;   % seconds

figure;
hold on 
plot(t*1e6, frame_1, 'LineWidth', 1.2);  % Time in microseconds
xlabel('Time (\mus)');
ylabel('Amplitude');
title('Raw RF Signal Middle Channel Frame 1');
grid on;
%}

%% Conversion to frequency domain 

%nfft = 2^18;
[pxx, f] = pwelch(frame_1, [], [], [], fs);

% average magnitude across channels
      
pxx_db = 10*log10(pxx)

% average magnitude across frames for one channel 
%mag_single_ch = mag2db(squeeze(X_single_ch)); 

% Plot 
plot(f/1e6, pxx_db, 'LineWidth', 1.5);  % Freq
hold on 
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
title('Probe Data Frequency Domain Middle Channel Frame 1');
xlim([0 14])  % 1.5x center frequency
ylim([-65, -40])  % Capture noise floor to peak
grid on;


%{
f0 = 0.5e6 % fundamental freq 

harmonics = 4:12; 

for h = harmonics
    h_freq = h * f0;
    [~, idx] = min(abs(f - h_freq));    % closest frequency bin
    plot(h_freq/1e6, pxx_db(idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    text(h_freq/1e6, pxx_db(idx)+5, sprintf('H%d', h), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

ultraharmonics = 4:12; 

for u = ultraharmonics
    u_freq = (u+0.5) * f0;
    [~, idx] = min(abs(f - u_freq));    % closest frequency bin
    plot(h_freq/1e6, pxx_db(idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    text(h_freq/1e6, pxx_db(idx)+5, sprintf('U%d', u), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

%}
%{
%% Extract for all datafiles 

% Process each datafile
for k = 1:length(data_files)
    file_path = string(fullfile(data_files(k).folder, data_files(k).name));
    data = imaging.common.acquisition_reader(file_path);
    disp(data.data_params)
end



% Save metadata
    all_data(k).filename = data_files(k).name;
    all_data(k).title = data.acq_title;
    all_data(k).timestamp = data.acq_timestamp;
    all_data(k).n_frames = data.number_of_frames;
    all_data(k).n_samples = data.rfdata_dimensions(1);
    all_data(k).n_channels = data.rfdata_dimensions(2);
    all_data(k).elapsed_time = data.elapsed_timevector;
% Save data in 3D array
    all_data(k).rf_data = zeros(data.rfdata_dimensions(1), ...
                                  data.rfdata_dimensions(2), ...
                                  data.number_of_frames);
    disp(all_data)
% Process each frame per file
    %for frame_idx = 1:data.number_of_frames
% Convert data from uint16 to double
        %all_data(k).rf_data(:, :, frame_idx) = double(data.data{frame_idx});
   % end
end
%save(fullfile(probe_path, 'processed_dataset2.mat'), 'all_data', '-v7.3');
%fprintf('All data saved\n');


file_idx = 1;  % First file
channel = 64;  % Middle channel
frame = 1;     % First frame


% Get sampling frequency and create time vector
fs = all_data(file_idx).fs;  % Sampling frequency in Hz
n_samples = all_data(file_idx).n_samples;
t = (0:n_samples-1) / fs;  % Time in seconds

fprintf('=== Plotting Data ===\n');
fprintf('File: %s\n', all_data(file_idx).filename);
fprintf('Title: %s\n', all_data(file_idx).title);
fprintf('Sampling frequency: %.2f MHz\n', fs/1e6);
fprintf('Channel: %d, Frame: %d\n\n', channel, frame);


rf_signal = all_data(file_idx).rf_data(:, channel, frame);

% Time vector for one frame
n_samples = acq_data.rfdata_dimensions(1);
t = (0:n_samples-1) / fs;  % Time in seconds

% Select middle channel and first frame
channel = 64;
frame = 1;
rf_signal = double(acq_data.data{frame}(:, channel));

% PLOT 1: Time domain signal (zoomed to see RF oscillations)
figure('Position', [100 100 1400 800]);

subplot(3,2,1);
plot(t*1e6, rf_signal);  % Time in microseconds
xlabel('Time (μs)');
ylabel('Amplitude');
title(sprintf('Full RF Signal - Channel %d, Frame %d', channel, frame));
grid on;

% Zoom into first 10 microseconds to see RF frequency
subplot(3,2,2);
zoom_samples = round(10e-6 * fs);  % First 10 microseconds
plot(t(1:zoom_samples)*1e6, rf_signal(1:zoom_samples), 'LineWidth', 1.5);
xlabel('Time (μs)');
ylabel('Amplitude');
title('Zoomed: First 10 μs (See RF Carrier Frequency)');
grid on;




% Example plots 
figure('Name', sprintf('File %d: %s - Frame 1', k, data.acq_title));
            
            % Plot 1: RF data as an image
            subplot(2,2,1);
            imagesc(rf_data_double);
            colorbar;
            title(sprintf('RF Data - Frame 1\n%s', data.acq_title));
            xlabel('Channel');
            ylabel('Sample (time)');
            
            % Plot 2: Single channel signal
            subplot(2,2,2);
            plot(rf_data_double(:, 64)); % Plot channel 64 (middle channel)
            title('Single Channel (64) RF Signal');
            xlabel('Sample');
            ylabel('Amplitude');
            grid on;

%}
