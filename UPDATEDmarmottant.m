clear; clc; close all;

% marmottant model function for a single bubble
% define simulation parameters
fs = 500e6;                 % sampling frequency [Hz]
T = 2e-6;                   % total simulation time [s]
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

% calculate PCD pressure signal
P = params.rho * (R .* Rddot + 2 * Rdot.^2); % far-field approx.

% plot time-domain PCD pressure
figure(1);
plot(t * 1e6, P);
xlim([0, T * 1e6]);
xlabel('Time [\mus]');
ylabel('Pressure [a.u.]');
title('PCD Signal from Single Bubble in Time Domain');

% FFT-based frequency domain analysis
N = length(P);
P_fft = fft(P);
P_fft_mag = abs(P_fft / N);
P_fft_mag = P_fft_mag(1:N/2+1);       % single-sided spectrum
P_fft_mag(2:end-1) = 2*P_fft_mag(2:end-1);

f = fs * (0:(N/2)) / N;               % frequency axis

% plot frequency-domain signal
figure(2);
plot(f * 1e-6, P_fft_mag); % convert to MHz
xlabel('Frequency [MHz]');
ylabel('Amplitude [a.u.]');
title('FFT of PCD Signal');

% plot radius oscillation
figure(3);
plot(t * 1e6, R * 1e6);
xlabel('Time [\mus]', 'FontSize', 12);
ylabel('Radius [\mum]', 'FontSize', 12);
title('Radius Oscillation (Marmottant Model)', 'FontSize', 14);

function params = setup_marmottant_params(R0, chi, A)
    params.R0 = R0;
    params.Rb = R0;
    params.chi = chi;
    params.sigma_break = 1.0;
    params.sigma_water = 0.072;
    params.rho = 1000;
    params.mu = 0.001;
    params.kappa_s = 7.2e-9;
    params.c = 1480;
    params.kappa = 1.095;
    params.P0 = 1e5;

    f = 2.9e6;
    cycles = 8;
    params.Pac = @(t) A * sin(2*pi*f*t) .* (t < cycles / f);
end

function dYdt = marmottant_ode(t, Y, p)
    R = Y(1); Rdot = Y(2); Pac = p.Pac(t);
    sigma = marmottant_surface_tension(R, p);
    sigma_R0 = marmottant_surface_tension(p.R0, p);
    Pg = (p.P0 + 2*sigma_R0/p.R0)*(p.R0/R)^(3*p.kappa);

    numerator = Pg * (1 - (3 * p.kappa * Rdot) / p.c) ...
              - p.P0 - 2*sigma/R ...
              - 4*p.mu*Rdot/R ...
              - 4*p.kappa_s*Rdot/R^2 ...
              - Pac;
    denominator = p.rho * R;

    Rddot = numerator / denominator - (3/2) * Rdot^2 / R;
    dYdt = [Rdot; Rddot];
end

function sigma = marmottant_surface_tension(R, p)
    Rb = p.Rb;
    Rrupture = Rb * sqrt(1 + p.sigma_break / p.chi);
    if R <= Rb
        sigma = 0;
    elseif R <= Rrupture
        sigma = 2 * p.chi * (R / Rb - 1);
    else
        sigma = p.sigma_water;
    end
end
