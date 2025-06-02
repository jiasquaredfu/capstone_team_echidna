clear, clc

% marmottant model function for a single bubble
% define simulation parameters
fs = 50e6;                 % sampling frequency [Hz]
T = 2e-6;                  % total simulation time [s]
t = linspace(0, T, fs*T);  % time vector

% define bubble parameters
R0 = 0.975e-6;             % equilibrium radius [m]
chi = 1.0;                 % shell elasticity [N/m]
A = 1.3e5;                 % acoustic pressure amplitude [Pa]

% call function to setup parameters
params = setup_marmottant_params(R0, chi, A);

% define initial conditions and solve ODE
Y0 = [R0, 0]; 
[~, Y] = ode45(@(t, Y) marmottant_ode(t, Y, params), t, Y0);


% extract bubble dynamics
R = Y(:, 1);

fprintf('Min R: %.6f µm\\n', min(R)*1e6);
fprintf('Max R: %.6f µm\\n', max(R)*1e6);
fprintf('Amplitude: %.6f µm\\n', (max(R) - min(R))*1e6);

Rdot = Y(:, 2);
Rddot = gradient(Rdot, t);

% calculate PCD pressure signal
P = params.rho * (R .* Rddot + 2 * Rdot.^2);

% add small Gaussian noise
P = P + 0.05 * max(abs(P)) * randn(size(P));

% plot PCD pressure signal in time domain
figure;
plot(t * 1e6, P);
xlim([0, T * 1e6]);
xlabel('Time [\mus]');
ylabel('Pressure [a.u.]');
title('PCD Signal from Single Bubble in Time Domain');

% use pwelch to apply FFT and plot in frequency-domain
figure;
pwelch(P, [], [], [], fs, 'power');
title('FFT PCD Signal in Frequency Domain');
xlabel('Frequency [Hz]');

% try to replicate figure 5b from paper
figure;
plot(t * 1e6, R * 1e6, 'k', 'LineWidth', 1.25);
xlabel('Time [\\mus]', 'FontSize', 12);
ylabel('Radius [\\mum]', 'FontSize', 12);
title('Radius Oscillation (Marmottant Model)', 'FontSize', 14);
xlim([0, T * 1e6]); ylim([0.8, 1.05]);
grid on; set(gca, 'FontSize', 12);


% supporting functions
% function to setup marmottant model parameters
function params = setup_marmottant_params(R0, chi, A)
    params.R0 = R0;                   % equilibrium radius
    params.Rb = R0;                   % buckling radius
    params.chi = chi;                 % shell elasticity modulus
    params.sigma_break = 1.0;         % surface tension rupture threshold
    params.sigma_water = 0.072;       % surface tension of water
    params.rho = 1000;                % liquid density
    params.mu = 0.001;                % liquid viscosity
    params.mu_s = 15e-9;              % shell surface viscosity/dilatational viscosity 
    params.c = 1480;                  % speed of sound in liquid
    params.gamma = 1.095;             % polytropic exponent
    params.P0 = 1e5;                  % ambient pressure 

    f = 3e6;                          % frequency
    cycles = 8;                       % number of cycles
    params.Pac = @(t) A * sin(2*pi*f*t) .* (t < cycles / f);
end

% define marmottant ODE function
function dYdt = marmottant_ode(t, Y, p)
    R = Y(1);
    Rdot = Y(2);
    Pac = p.Pac(t);
    Pg = p.P0 * (p.R0 / R)^(3 * p.gamma) * (1 - 3 * Rdot / p.c);
    sigma = marmottant_surface_tension(R, p);

    Rddot = (Pg - p.P0 - 2 * sigma / R ...
             - 4 * p.mu * Rdot / R ...
             - 4 * p.mu_s * Rdot / R^2 ...
             - Pac) / p.rho ...
             - (3/2) * Rdot^2 / R;

    dYdt = [Rdot; Rddot];
end

% define marmottant surface tension function
function sigma = marmottant_surface_tension(R, p)
    Rb = p.Rb;
    Rrupture = Rb * sqrt(1 + p.sigma_break / p.chi);

    if R <= Rb
        sigma = 0;
    elseif R <= Rrupture
        sigma = p.chi * ((R^2 / Rb^2) - 1);
    else
        sigma = p.sigma_water;
    end
end