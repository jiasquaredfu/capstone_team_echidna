%% Dataset 2 Analysis 
% Capstone II - Team E(Chidna)
% JiaJia Fu 
% 3/14/26

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

figure;
grid on
hold on

f0 = 0.5e6;
f_norm = f / f0;

% ---- Full Spectrum ----
h_full = plot(f_norm, pxx_db, 'Color',[0.5, 0.5, 0.5], 'LineWidth', 1.5, ...
              'DisplayName','Full PSD');


% ---- Harmonics ----
[hf, hv, uf, uv] = get_harmonics(f, pxx_db);

hf_norm = hf / f0;
uf_norm = uf / f0;

h_harm = plot(hf_norm, hv, 'ro', ...
              'MarkerFaceColor','r', ...
              'MarkerSize',8, ...
              'DisplayName','Harmonics');

% ---- Label Harmonics ----
for i = 1:length(hf_norm)
    text(hf_norm(i), hv(i)+1, sprintf('H%d',i+5), ...
        'HorizontalAlignment','center', ...
        'FontWeight','bold');
end

% ---- Ultraharmonics ----
h_ultra = plot(uf_norm, uv, '^', ...
               'Color',[0 0.5 0], ...
               'MarkerFaceColor',[0 0.5 0], ...
               'MarkerSize',8, ...
               'DisplayName','Ultraharmonics');

for i = 1:length(uf_norm)
    text(uf_norm(i), uv(i)+1, sprintf('U%.1f',6.5+i-1), ...
        'HorizontalAlignment','center', ...
        'FontWeight','bold');
end

% ---- Broadband Noise Mask ----
[noise_mask, noise_psd_db, broadband_noise_db] = get_bbnoise(f, pxx, hf, uf)
pxx_noise_db = pxx_db;
pxx_noise_db(~noise_mask) = NaN;

h_noise = plot(f_norm, pxx_noise_db, 'b','LineWidth',2);

% ---- Baseline ----
baseline = median(noise_psd_db); % more accurate than mean?

h_base = yline(baseline, '--', ...
               'Color',[0, 0, 0], ...
               'LineWidth',2, ...
               'DisplayName','Baseline');

% ---- Labels ----
xlabel('Normalized Frequency (f/f_0)', ...
       'FontSize',20,'FontWeight','bold')

ylabel('Magnitude (dB)', ...
       'FontSize',20,'FontWeight','bold')

title('Probe Data Frequency Domain - Middle Channel Frame 1', ...
      'FontSize',20)

xlim([5 13])
ylim([-65 -40])

ax = gca;
ax.FontSize = 14;

% ---- Legend ----
legend([h_full h_noise h_base h_harm h_ultra], ...
       {'Full PSD','Broadband Region','Baseline', ...
        'Harmonics','Ultraharmonics'}, ...
       'FontSize',14, ...
       'Location','best')

hold off



%% Extract for all datafiles

data_files = dir(fullfile(probe_path, '*.pacq*'));
f0         = 0.5e6;
ALL_CHANNELS = 128;

% --- Preallocation ---
temp              = extract_probe_data(probe_path, 1);
N_frames_per_file = size(temp.rf_data, 3);
total_rows        = length(data_files) * N_frames_per_file;

% Build dummy feature struct to type the preallocated array
dummy_f   = linspace(0, 10e6, 1000)';
dummy_pxx = ones(size(dummy_f));
dummy_feat = get_features(dummy_f, dummy_pxx, f0);  % FIX: was build_feature_row (deleted)
dummy_feat.filename = 'dummy.pacq';                 % FIX: add metadata fields to dummy
dummy_feat.frame    = 0;

all_features(total_rows, 1) = dummy_feat;
row_idx = 1;

% Main loop
for i = 1:length(data_files)

    clean_data = extract_probe_data(probe_path, i);  % FIX: removed stray semicolon-less call
    N_frames   = size(clean_data.rf_data, 3);

    for frame_idx = 1:N_frames

        pxx_sum = [];

        for ch = 1:ALL_CHANNELS
            frame = clean_data.rf_data(:, ch, frame_idx);
            [pxx_ch, f] = pwelch(frame, [], [], [], clean_data.fs);

            if isempty(pxx_sum)
                pxx_sum = zeros(size(pxx_ch));
            end
            pxx_sum = pxx_sum + pxx_ch;
        end

        pxx_avg = pxx_sum / ALL_CHANNELS;  % average in linear power domain

        feat = get_features(f, pxx_avg, f0);

        feat.filename = data_files(i).name;
        feat.frame    = frame_idx;

        all_features(row_idx) = feat;
        row_idx = row_idx + 1;
    end
end

% Trim in case N_frames varied across files
all_features = all_features(1:row_idx-1);

% Convert and save
features_table = struct2table(all_features);
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


function features = get_features(f, pxx, f0)
    pxx_db = 10*log10(pxx);
    df     = f(2) - f(1);

    harmonics      = 6:12;
    ultraharmonics = 6.5:1:12.5;
    search_window  = round(0.08e6 / df);

    features = struct();

    % Accumulators for BB noise exclusion — populated during loops below
    harm_freq_vals  = zeros(1, length(harmonics));
    ultra_freq_vals = zeros(1, length(ultraharmonics));

    % --- Harmonics ---
    for i = 1:length(harmonics)
        target_f  = harmonics(i) * f0;
        [~, cidx] = min(abs(f - target_f));
        idx_range = max(1, cidx-search_window) : min(length(f), cidx+search_window);

        [~, local] = max(pxx_db(idx_range));
        peak_idx   = idx_range(local);

        h = harmonics(i);
        features.(sprintf('H%d_freq',  h)) = f(peak_idx);
        features.(sprintf('H%d_db',    h)) = pxx_db(peak_idx);
        features.(sprintf('H%d_power', h)) = trapz(f(idx_range), pxx(idx_range));

        harm_freq_vals(i) = f(peak_idx);   % save for BB exclusion
    end

    % --- Ultraharmonics ---
    for i = 1:length(ultraharmonics)
        target_f  = ultraharmonics(i) * f0;
        [~, cidx] = min(abs(f - target_f));
        idx_range = max(1, cidx-search_window) : min(length(f), cidx+search_window);

        [~, local] = max(pxx_db(idx_range));
        peak_idx   = idx_range(local);

        u = ultraharmonics(i) * 10;
        features.(sprintf('U%d_freq',  u)) = f(peak_idx);
        features.(sprintf('U%d_db',    u)) = pxx_db(peak_idx);
        features.(sprintf('U%d_power', u)) = trapz(f(idx_range), pxx(idx_range));

        ultra_freq_vals(i) = f(peak_idx);  % save for BB exclusion
    end

    % --- Broadband noise ---
    all_peak_freqs = [harm_freq_vals(:); ultra_freq_vals(:)];
    noise_range_min = 9 * f0;
    noise_range_max = 10.5 * f0;
    noise_exclusion_Hz = 0.08e6;

    noise_mask = (f >= noise_range_min) & (f <= noise_range_max);
    for k = 1:length(all_peak_freqs)
        noise_mask(abs(f - all_peak_freqs(k)) <= noise_exclusion_Hz) = false;
    end

    noise_psd = pxx(noise_mask);
    features.broadband_noise_db = 10*log10(mean(noise_psd));
    features.broadband_noise_power = trapz(f(noise_mask), noise_psd);
end


