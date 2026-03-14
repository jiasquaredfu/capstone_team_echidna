%% Dataset 2 Analysis 
% Capstone II - Team E(Chidna)
% JiaJia Fu 
% 2/17/26

clear;
clc;
close all;


% Add path of imaging function 
function_path = '/Users/jiasquared/Desktop/CODING/capstone_team_echidna/04_Projects/';
addpath(function_path);

% Add path of Philips-ATL L12-5 Probe Data to read 
probe_path = ['/Users/jiasquared/Desktop/CODING/capstone_team_echidna' ...
    '/1_and_2_2026-01-16_single_tube_sweep/probe'];


clean_data = extract_probe_data(probe_path, 1)
fs = clean_data.fs


%% Conversion to frequency domain for middle channel frame 1 example 

frame_1 = clean_data.rf_data(:, 64, 1); % middle channels frame 1 sample 
[pxx, f] = pwelch(frame_1, [], [], [], fs);

% convert to decibels for display
pxx_db = 10*log10(pxx)


%% Plot for middle channel frame 1 example 

f0 = 0.5e6; 
window_f = 0.3e6;
window_uf = 0.05e6;  
harmonics = 6:20; 

% Normalized Data Graph 
figure;
f_norm = f/f0
plot(f_norm, pxx_db, 'LineWidth', 1.5);  % Freq
hold on 
xlabel('Normalized Frequency');
ylabel('Magnitude (dB)');
title('Probe Data Frequency Domain Middle Channel Frame 1');
xlim([0 20])  % 1.5x center frequency
ylim([-65, -40])  % Capture noise floor to peak


% Data with Harmonics and Ultraharmonics
figure;
plot(f/1e6, pxx_db, 'LineWidth', 1.5);  % Freq
hold on 
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');
title('Probe Data Frequency Domain Middle Channel Frame 1');
xlim([0 10])  % 1.5x center frequency
ylim([-65, -40])  % Capture noise floor to peak


for h = harmonics
    target_f = h * f0;
    h_idx_range = find(f >= target_f-window_f & f <= target_f+window_f);

    [max_val, id_local] = max(pxx(h_idx_range));
    peak_idx = h_idx_range(id_local);

    plot(f(peak_idx)/1e6, pxx_db(peak_idx), ...
         'ro', 'MarkerFaceColor','r')
     text(f(peak_idx)/1e6, pxx_db(peak_idx)+3, ...
         sprintf('H%d', h), ...
         'HorizontalAlignment','center', ...
         'FontWeight','bold');

end


for h = harmonics
    target_f = (h+0.5) * f0;
    u_idx_range = find(f >= target_f-window_uf & f <= target_f+window_uf);

    [max_val, id_local] = max(pxx(u_idx_range));
    peak_idx = u_idx_range(id_local);


    plot(f(peak_idx)/1e6, pxx_db(peak_idx), ...
         'go', 'MarkerFaceColor','g')
     text(f(peak_idx)/1e6, pxx_db(peak_idx)+3, ...
         sprintf('U%d', h), ...
         'HorizontalAlignment','center', ...
         'FontWeight','bold');

end


%% Get broadband noise 

noise_exclusion_Hz = 0.05e6;  

% Create mask for all frequency bins
noise_mask = true(size(pxx));

% Combine all harmonic and ultraharmonic target frequencies
h_freqs = harmonics * f0;
u_freqs = (harmonics + 0.5) * f0;
all_peak_freqs = [h_freqs, u_freqs];

% Exclude regions around each peak using frequency windows
for k = 1:length(all_peak_freqs)
    freq_center = all_peak_freqs(k);
    noise_mask(f >= (freq_center - noise_exclusion_Hz) & f <= (freq_center + noise_exclusion_Hz)) = false;
end

% Extract PSD values in noise regions (mV^2/Hz)
noise_psd_values = pxx(noise_mask);

% Broadband noise (mean PSD)
broadband_noise_mean_psd = mean(noise_psd_values);
broadband_noise_mean_psd_dB = 10*log10(broadband_noise_mean_psd);


figure;
set(gcf,'Position',[100 100 1800 1000])
hold on;

% Full spectrum (broadband)
h_full = plot(f/1e6, 10*log10(pxx), 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'DisplayName','Broadband Spectrum');

% Harmonics
for h = harmonics
    target_f = h * f0;
    h_idx_range = find(f >= target_f-window_f & f <= target_f+window_f);
    [~, id_local] = max(pxx(h_idx_range));
    peak_idx = h_idx_range(id_local);

    plot(f(peak_idx)/1e6, pxx_db(peak_idx), 'ro', 'MarkerFaceColor','r', 'MarkerSize',15);  % no DisplayName
    text(f(peak_idx)/1e6, pxx_db(peak_idx)+3, sprintf('H%d', h), ...
        'color','red', 'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize', 18);
end

% Ultraharmonics
for h = harmonics
    target_f = (h+0.5) * f0;
    u_idx_range = find(f >= target_f-window_uf & f <= target_f+window_uf);
    [~, id_local] = max(pxx(u_idx_range));
    peak_idx = u_idx_range(id_local);

    plot(f(peak_idx)/1e6, pxx_db(peak_idx), 'go','color',[0 0.5 0], 'MarkerFaceColor','[0 0.5 0]', 'MarkerSize',15);  % no DisplayName
    text(f(peak_idx)/1e6, pxx_db(peak_idx)+3, sprintf('U%d', h), ...
        'color',[0 0.5 0],'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize', 18);
end

% Noise regions
h_noise = plot(f(noise_mask)/1e6, 10*log10(pxx(noise_mask)), 'b','DisplayName','Noise Regions');

% Labels and limits
xlabel('Normalized Frequency (f/f0)', 'FontSize', 30, 'FontWeight','bold');
ylabel('Magnitude (dB)', 'FontSize', 30, 'FontWeight','bold');
title('Probe Data Frequency Domain Middle Channel Frame 1', 'FontSize', 30);
ax = gca; % Get current axes
ax.FontSize = 16; 
xlim([0 8]);
ylim([-65 -40]);

% Legend with only the desired two entries
legend([h_full, h_noise], 'FontSize', 20, 'FontWeight','bold');
%% Extract for all datafiles 

data_files = dir(fullfile(probe_path, '*.pacq*'));

all_features = [];  % struct array to store features
channels = [1, 64, 128];  % first, middle, last channels

% Process each data file
for k = 1:length(data_files)
    
    % Extract RF data
    clean_data = extract_probe_data(probe_path, k);
    N_frames = size(clean_data.rf_data, 3);
    
    % Loop over selected channels
    for ch = channels
        % Loop over all frames
        for frame_idx = 1:N_frames
            frame = clean_data.rf_data(:, ch, frame_idx);  % single column vector
            
            % Compute PSD
            [pxx, f] = pwelch(frame, [], [], [], clean_data.fs);
            pxx_db = 10*log10(pxx);
            
            % Extract harmonic features
            features = get_features(f, pxx, pxx_db, f0);
            
            
            % Add metadata
            features.filename = data_files(k).name;
            features.channel = ch;
            features.frame = frame_idx;
            
            % Append to results
            all_features = [all_features; features];
        end
    end
end

% Convert struct array to table
features_table = struct2table(all_features);

% Save to CSV
writetable(features_table, 'probe_data_all_frames.csv');

disp('Feature extraction complete for all frames!');



%% Functions 

function out = extract_probe_data(probe_path, i)
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
        file_path = string(fullfile(data_files(i).folder, data_files(i).name));
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


function features = get_features(f, pxx, pxx_db, f0)
    noise_exclusion_Hz = 0.05e6;  
    harmonics = 6:20;
    window_f = 0.3e6;
    window_uf = 0.05e6;

    features = struct();

    % Harmonics
    for h = harmonics
        target_f = h * f0;
        idx_range = find(f >= target_f-window_f & f <= target_f+window_f);

        if isempty(idx_range), continue; end

        [~, id_local] = max(pxx(idx_range));
        peak_idx = idx_range(id_local);

        % calculate power from psd by taking integral 
        features.(['H',num2str(h),'_freq']) = f(peak_idx);
        features.(['H',num2str(h),'_power']) = trapz(f(idx_range), pxx(idx_range));
        features.(['H',num2str(h),'_db']) = pxx_db(peak_idx);
    end

    % Ultra-harmonics
    for h = harmonics
        target_f = (h + 0.5) * f0;
        idx_range = find(f >= target_f-window_uf & f <= target_f+window_uf);

        if isempty(idx_range), continue; end

        [~, id_local] = max(pxx(idx_range));
        peak_idx = idx_range(id_local);

        features.(['U',num2str(h),'_freq']) = f(peak_idx);
        features.(['U',num2str(h),'_power']) = trapz(f(idx_range), pxx(idx_range));
        features.(['U',num2str(h),'_db']) = pxx_db(peak_idx);
    end

        %Broad-band Noise 
    
        % Create mask for all frequency bins
    noise_mask = true(size(pxx));

    noise_range_min = 1.5;  
    noise_range_max = 3.0;  

    
    % Combine all harmonic and ultraharmonic target frequencies
    h_freqs = harmonics * f0;
    u_freqs = (harmonics + 0.5) * f0;
    all_peak_freqs = [h_freqs, u_freqs];
    
    % Exclude regions around each peak using frequency windows
    for k = 1:length(all_peak_freqs)
        freq_center = all_peak_freqs(k);
        noise_mask(f >= (freq_center - noise_exclusion_Hz) & f <= (freq_center + noise_exclusion_Hz)) = false;
    end
    
    % Extract PSD values in noise regions (mV^2/Hz)
    noise_psd_values = pxx(noise_mask);
    
    % Store mean PSD and dB
    features.broadband_noise_psd = mean(noise_psd_values);
    features.broadband_noise_dB  = 10*log10(features.broadband_noise_psd);

end
