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

N = length(signal);           % number of samples
y = fft(signal);              % compute FFT

% --- Frequency vector (positive only) ---
f = (0:N/2)*(fs/N);          

% --- Normalize and compute one-sided amplitude spectrum ---
P2 = abs(y/N);                % full amplitude spectrum
P1 = P2(1:N/2+1);             % take only positive half
P1(2:end-1) = 2*P1(2:end-1);

% --- Plot ---
figure;
plot(f,P1)
title('Vokurka fft graph')
xlabel("Frequency (Hz)")
ylabel("Amplitude [PA]")

%%

n_samples = 5000;      % Number of samples
f0 = 1e6;              % Pure tone frequency (1 MHz)
A = 1;                 % Amplitude

% Time vector
t = (0:n_samples-1)/fs;

% Pure sine wave at 1 MHz
signal = A * sin(2*pi*f0*t);

N = length(signal);
y = fft(signal);
f = (0:N/2)*(fs/N);             % Frequency vector (positive side)
P2 = abs(y/N);                  % Two-sided amplitude spectrum
P1 = P2(1:N/2+1);               % One-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);    % Energy correction

% Plot
figure("Name", "FFT of Pure 1 MHz Tone", "Units", "inches", "Position", [4 4 8 4]);
plot(f, P1);
xlim([0 5e6]);  % Zoom in to see around 1 MHz
xlabel("Frequency [Hz]");
ylabel("|Y(f)|");
title("FFT of Pure 1 MHz Sine Wave");
