% clear, clc
% modified marmottant model
%% main simulation funciton
% define function to generate pcd signal using marmottant model
function marmottant
    % define simulation parameters
    N_bubbles = 50;            % number of individual bubbles
    fs = 50e6;                 % sampling frequency [Hz]
    T = 5e-6;                  % total simulation duration [s]
    t = linspace(0, T, fs*T);  % time vector

    % initialize empty signal for accumulating pressure from bubbles
    P_total = zeros(size(t)); 

    % loop over each bubble
    for i = 1:N_bubbles
        % randomize parameters
        R0 = normrnd(2e-6, 0.2e-6);                          % equilibrium radius where mean=2 [m]
        chi = normrnd(1.0, 0.2);                             % shell elasticity [N/m]
        A = normrnd(1e5, 2e4);                               % acoustic pressure amplitude [Pa]
        shift_samples = randi([0, round(0.1 * length(t))]);  % random phase/time shift

        % setup bubble parameters
        params = setup_marmottant_params(R0, chi, A);
        Y0 = [R0, 0]; % initial condition=bubble at rest

        % solve ode to simulate marmottant dynamics
        [~, Y] = ode45(@(t, Y) marmottant_ode(t, Y, params), t, Y0);

        % extract radius
        R = Y(:,1);

        % compute acceleration using gradient
        Rdot = Y(:,2);
        Rddot = gradient(Rdot, t);

        % compute pressure signal emitted at bubble wall
        P = params.rho * (R .* Rddot + 2 * Rdot.^2);

        % random amplitude scaling and time shift
        P = circshift(P, shift_samples);
        P = P * rand();

        % add bubble signal to total PCD signal
        P_total = P_total + P(:)';
    end

    % add 5% white Gaussian noise 
    % this simulates realistic detector noise
    noise_level = 0.05 * max(abs(P_total));
    P_total = P_total + noise_level * randn(size(P_total));

    % plot PCD signal
    figure;
    plot(t * 1e6, P_total);
    xlabel('Time [\mus]');
    ylabel('Pressure [a.u.]');
    title('Synthetic PCD-like Signal from Modified Marmottant Model');
end
%% define funtion to setup marmottant parameters
% returns shell & fluid properties, acoustic pulse definition
function params = setup_marmottant_params(R0, chi, A)
    params.R0 = R0;
    params.Rb = 0.975 * R0;
    params.chi = chi;
    params.sigma_break = 0.1;
    params.sigma_water = 0.072;
    params.rho = 1000;
    params.mu = 0.001;
    params.mu_s = 1e-8;
    params.c = 1480;
    params.gamma = 1.4;
    params.P0 = 1e5;

    f = 3e6; cycles = 8;
    params.Pac = @(t) A * sin(2 * pi * f * t) .* (t < cycles / f);
end
%% define differential equations
% use Rayleigh-Plesset model with surface tension from Marmottant model
% mu=liquid viscosity, mu_s=shell surface viscosity
function dYdt = marmottant_ode(t, Y, p)
    R = Y(1);
    Rdot = Y(2);
    Pac = p.Pac(t);
    Pg = p.P0 * (p.R0 / R)^(3 * p.gamma); 
    sigma = smooth_surface_tension(R, p);
    damping_radiation = (4/3) * pi * R^2 * p.rho / p.c; % add radiation damping

    Rddot = (Pg - p.P0 - Pac ...
             - (2 * sigma / R) ...
             - (4 * p.mu * Rdot / R) ...
             - (4 * p.mu_s * Rdot / R^2)) ...
             / (p.rho + damping_radiation) ...
             - (3/2) * (Rdot^2 / R);

    dYdt = [Rdot; Rddot];
end
%% define function to compute surface tension
% considers buckling, elastic region, and rupture states
% smooth transitions using sigmoid curves
function sigma = smooth_surface_tension(R, p)
    delta = 0.01 * p.Rb;
    sigma_elastic = p.chi * (R.^2 / p.Rb^2 - 1);
    w_buckled = 1 ./ (1 + exp((R - p.Rb) / delta));
    Rrupture = p.Rb * sqrt(1 + p.sigma_break / p.chi);
    w_rupture = 1 ./ (1 + exp((Rrupture - R) / delta));

    sigma = w_buckled .* 0 + ...
            (1 - w_buckled - w_rupture) .* sigma_elastic + ...
            w_rupture .* p.sigma_water;
end