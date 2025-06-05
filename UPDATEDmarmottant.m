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


% try to replicate figure 5b from paper
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
Fs = 500e6;
%Compute the Fourier transform of the signal. Plot the magnitude of the transform as a function of frequency.

y = fft(R);
ly = length(y);
f = (-ly/2:ly/2-1)/ly*Fs;
figure(2)
plot(f,abs(y))
xlim([0 2.5e8])
title("Double-Sided Amplitude Spectrum of x(t)")
xlabel("Frequency (Hz)")
ylabel("|y|")
grid