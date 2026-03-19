clear, clc, close all;

% marmottant model function for a single bubble
% define simulation parameters
fs = 500e6;                 % sampling frequency [Hz]
T = 2e-6;                  % total simulation time [s]
t = linspace(0, T, fs*T);  % time vector 

% define bubble parameters
R0 = 0.975e-6;             % equilibrium radius [m]
chi = 0.55;                % shell elasticity [N/m]
A = 130e3;                 % acoustic pressure amplitude [Pa]

% call function to setup parameters
params = setup_marmottant_params(R0, chi, A);

% define initial conditions and solve ODE
Y0 = [R0, 0]; 
[~, Y] = ode45(@(t, Y) marmottant_ode(t, Y, params), t, Y0);

% extract bubble dynamics
R = Y(:, 1);
Rdot = Y(:, 2); % [m/s]
Rddot = gradient(Rdot, t); % [m/s^2]

% calculate PCD pressure signa
P = params.rho * (R .* Rddot + 2 * Rdot.^2); % pressure detected by PCD far-field approx.

% plot PCD pressure signal in time domain
figure(1);
plot(t * 1e6, P, 'o-');
xlim([0, T * 1e6]);
xlabel('Time [\mus]');
ylabel('Pressure [Pa]');
title('PCD Signal from Single Bubble in Time Domain');

% FFT of PCD signal
L = length(P);              % length of signal
Yf = fft(P);                % compute FFT
f = fs * (0:(L/2)) / L;     % (positive) frequency vector

% get amplitude spectrum
P_amp = abs(Yf / L); % divide FFT by L to normalize amplitude
P_amp = P_amp(1:L/2+1); % take positive frequencies only
% since we only take the first half of the frequencies
P_amp(2:end-1) = 2 * P_amp(2:end-1);  % double remaining amplitudes except for 0Hz and fs/2

% plot fft
figure(2);
% plot amplitude vs frequency
plot(f / 1e6, P_amp, 'LineWidth', 1.25);  % convert frequency to MHz
xlabel('Frequency [MHz]');
ylabel('Amplitude [a.u.]');
title('FFT of PCD Signal');
grid on;
xlim([0 50]);         

% replicate figure 5b from paper
figure(3);
plot(t * 1e6, R * 1e6, 'k', 'LineWidth', 1.25, color='blue');
xlabel('Time [\mus]', 'FontSize', 12);
ylabel('Radius [\mum]', 'FontSize', 12);
title('Radius Oscillation (Marmottant Model)', 'FontSize', 14);

% supporting functions
% function to setup marmottant model parameters
function params = setup_marmottant_params(R0, chi, A)
    params.R0 = R0;                   % equilibrium radius [m]
    params.Rb = R0;                   % buckling radius [m]
    params.chi = chi;                 % shell elasticity modulus [N/m]
    params.sigma_break = 1.0;         % surface tension rupture threshold [N/m]
    params.sigma_water = 0.072;       % surface tension of water [N/m]
    params.rho = 1000;                % liquid density [kg/m^3]
    params.mu = 0.001;                % liquid viscosity [Pa*s = N*s/m^2]
    params.kappa_s = 7.2e-9;          % surface dilatational viscosity [N/m]
    params.c = 1480;                  % speed of sound in liquid [m/s]
    params.kappa = 1.095;             % polytropic exponent [dimensionless]
    params.P0 = 1e5;                  % ambient pressure [Pa]

    f = 2.9e6;                        % frequency [Hz]
    cycles = 8;                       % number of cycles
    params.Pac = @(t) A * sin(2*pi*f*t) .* (t < cycles / f); % sinusoidal ultrasound acoustic pressure
end

% define marmottant ODE function
function dYdt = marmottant_ode(t, Y, p)
    R = Y(1);
    Rdot = Y(2);
    Pac = p.Pac(t);         % acoustic pressure at time t
    
    sigma = marmottant_surface_tension(R, p); % surface tension at current radius
    sigma_R0= marmottant_surface_tension(p.R0, p); % surface tension at equilibrium radius

    Pg = (p.P0 + 2 * sigma_R0 / p.R0) * (p.R0 / R)^(3 * p.kappa); % gas pressure from polytropic law

    numerator = Pg * (1 - (3 * p.kappa * Rdot) / p.c) ...
            - p.P0 ...
            - (2 * sigma) / R ...
            - (4 * p.mu * Rdot) / R ...
            - (4 * p.kappa_s * Rdot) / R^2 ...
            - Pac;

    denominator = p.rho * R;

    Rddot = numerator / denominator - (3/2) * Rdot^2 / R;
 
    dYdt = [Rdot; Rddot];
end

% define marmottant surface tension function
function sigma = marmottant_surface_tension(R, p)
    Rb = p.Rb;
    Rrupture = Rb * sqrt(1 + p.sigma_break / p.chi);

    if R <= Rb
        sigma = 0;
    elseif R <= Rrupture
        sigma = 2* p.chi * (R / Rb - 1);
    else
        sigma = p.sigma_water;
    end
end
%%
% manual windowing
f0 = 2.9e6;                      % fundamental frequency [Hz]
harmonics = [1 2 3] * f0;
ultraharmonics = [1.5 2.5] * f0;
window_half_width = 100e3;       % ±100 kHz window
threshold = 0.02 * max(P_amp);   % relative threshold for detection
labels = {};

fprintf('\nManual Frequency Window Analysis:\n');

% harmonics
for i = 1:length(harmonics)
    f_center = harmonics(i);
    idx = f >= (f_center - window_half_width) & f <= (f_center + window_half_width);
    amp = P_amp(idx);
    peak_val = max(amp);
    labels{end+1} = sprintf("Harmonic %d peak: %.4f", i, peak_val);
    fprintf("Harmonic %d: Peak = %.4f\n", i, peak_val);
end

% ultraharmonics
for i = 1:length(ultraharmonics)
    f_center = ultraharmonics(i);
    idx = f >= (f_center - window_half_width) & f <= (f_center + window_half_width);
    amp = P_amp(idx);
    peak_val = max(amp);
    labels{end+1} = sprintf("Ultraharmonic %.1f peak: %.4f", i + 0.5, peak_val);
    fprintf("Ultraharmonic %.1f: Peak = %.4f\n", i + 0.5, peak_val);
end

% cavitation classification
has_ultraharmonics = false;
has_harmonics = false;

for i = 1:length(labels)
    label = labels{i};

    if contains(label, 'Harmonic') || contains(label, 'Ultraharmonic')
        parts = strsplit(label, ':');
        if length(parts) == 2
            val = str2double(strtrim(parts{2}));
            if isnan(val)
                continue;
            end
            if contains(label, 'Ultraharmonic') && val > threshold
                has_ultraharmonics = true;
            elseif contains(label, 'Harmonic') && val > threshold
                has_harmonics = true;
            end
        end
    end
end

if has_ultraharmonics
    cavitation_type = 'Inertial Cavitation';
    box_color = [1, 0.85, 0.85]; % light red
elseif has_harmonics
    cavitation_type = 'Stable Cavitation';
    box_color = [0.85, 1, 0.85]; % light green
else
    cavitation_type = 'No Cavitation';
    box_color = [0.9, 0.9, 0.9]; % gray
end

fprintf("\nCavitation Type: %s\n", cavitation_type);

% plot FFT PCD signal with cavitation labels
figure(4);
plot(f / 1e6, P_amp, 'b', 'LineWidth', 1.5);
xlabel('Frequency [MHz]');
ylabel('Amplitude [a.u.]');
title('FFT of PCD Signal with Manual Window Labels');
grid on;
xlim([0 15]);

hold on;

% define window parameters
window_centers = [2.9, 5.8, 8.7, 4.35, 7.25];  % in MHz
window_half_width = 0.05;                      % in MHz
window_labels = {'Harmonic 1', 'Harmonic 2', 'Harmonic 3', ...
                 'Ultraharmonic 1.5', 'Ultraharmonic 2.5'};

% shade windows and add labels
for i = 1:length(window_centers)
    x1 = window_centers(i) - window_half_width;
    x2 = window_centers(i) + window_half_width;

    % shade window
    fill([x1 x2 x2 x1], [0 0 max(P_amp)*1.1 max(P_amp)*1.1], ...
         [0.85 0.85 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % add window text
    text(window_centers(i), max(P_amp)*1.05, window_labels{i}, ...
         'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end

% add cavitation label
annotation('textbox', [0.18, 0.78, 0.2, 0.06], ...  % smaller box
    'String', 'Cavitation Type: Stable Cavitation', ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'BackgroundColor', [0.9 1 0.9], ...
    'EdgeColor', 'k', ...
    'HorizontalAlignment', 'center');

hold off;