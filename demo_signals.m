            %% Example script on the Model system
            % N. Caso 2024
            clear; clc; close all;

            % MATLAB runs in the directory that you see in the panel to the
            % left. The functions are all in the folder called 03_Classes.
            % Let's add that to the PATH variable, which contains a list of
            % all the functions MATLAB knows:
            addpath("03_Classes\");
            
% 2 convenient variables for tracking units, MegaHertz and MegaPascal
MHz         = 1e6;
MPa         = 1e6;

% Set a sample rate. Needed to define the signals.
fs          = 100*MHz;

            %% Sinewave demo:

n_samples   = 500;
drive_freq  = 0.5*MHz;
amplitude   = 0.6*MPa;

            % this prints out the function documentation:
            help Sinusoid

            % MATLAB intellisense can help you with these functions too,
            % it will list out the arguments as you type.
sw          = Sinusoid(fs, n_samples, drive_freq, amplitude);
            
            % make a plot of the sinewave signal over time:
figure      (Name="Sinusoid Demo", Units="inches", Position=[1 1 8 3])
plot        (sw.time_vector, sw.time_series, "-k", LineWidth=1.1);
xlabel      ("Time [s]");
ylabel      ("Amplitude [PA]")

            %% Noise Demo:
noise       = Noise(fs, n_samples, 0.1*MPa, 1); % Generates random noise

            % make plots again:
sig         = sw + noise; % this creates a "MappedSignal" object, that contains the sum of the signals in its .time_series property.

hold        ("on");
% figure      (Name="Noise Demo", Units="inches", Position=[2 2 8 3])
plot        (sig.time_vector, sig.time_series, "--r", LineWidth=0.8);
xlabel      ("Time [s]");
ylabel      ("Amplitude [PA]")
legend      (["Clean Signal" "Noisy Signal"]);

            %% Vokurka Demo:
vok_signal  = Vokurka(fs, n_samples, ...
    10, ...         % number of sources
    drive_freq, ... % drive frequency of the FUS
    amplitude, ...  % source amplitude average
    0.1*MPa, ...    % source amplitude stanard deviation
    0, ...          % phase mean
    1e-7, ...       % phase standard deviation
    1e-7, ...       % time condition mean
    0.1e-9);        % time condition stanard deviation

            % plot
figure      (Name="Vokurka Demo", Units="inches", Position=[2 2 8 3])
plot        (vok_signal.time_vector, vok_signal.time_series, LineWidth= 1.1)
xlabel      ("Time [s]");
ylabel      ("Amplitude [PA]");

            %% Frequency Response Demo
            % an ultrasound probe doesn't have an infinite frequency
            % response. Let's say we're receiving with a Philips ATL L12-5.
            % make a filter that mimics this:
nyquist     = fs/2;
freq_limits = [5 12]*MHz./nyquist;

            % use a 4th order butterworth filter
[numerator, denominator] = butter(4, freq_limits, "bandpass");

            % make a plot of the filter frequency response (calling freqz
            % without assigning a variable does this)
figure; freqz(numerator, denominator); title("Filter Frequency Response")

            % plot again
figure      (Name="Filtered Vokurka Demo", Units="inches", Position=[3 3 8 3])
plot        (vok_signal.time_vector, vok_signal.filtfilt(numerator, denominator), LineWidth= 1.1)
xlabel      ("Time [s]");
ylabel      ("Amplitude [PA]");

%% Frequency-Domain Analysis using fft()

% --- Extract signal ---
signal = vok_signal.time_series;  

y = fft(signal);
N = length(signal);
f = (-N/2:N/2-1)/N*fs;
y_shifted = fftshift(y);
figure;
plot(f,abs(y_shifted))
title('Vokurka fft graph')
xlabel("Frequency (Hz)")
ylabel("Amplitude [PA]")

