% Run with EPG-X, revision 5da0e7128d2c7fc42ff13991ddeba806c81ed23f

% Sequence parameters
gamma = 267.5221 * 1e-3; % rad/ms/uT
B0 = 2.89 * 1e6; % uT
TR = 5; % ms
alpha = deg2rad(10); % rad
phi0 = 117; % deg

% Single-pool model
T1=779; % ms
T2=45; % ms

% MT model
T1_MT = [779 779]; % ms
T2_MT = 45; % ms
k_MT = 4.3e-3; % kHz
f_MT = 0.117; % unitless
G = 15.1; % us, absorption lineshape
B1 = 13; % uT

% Exchange model
T1x = [1000 500]; % ms
T2x = [100 20]; % ms
kx = 2e-3; % kHz
fx = 0.2; % unitless

delta_b = 2; % ppm
delta_b = delta_b*1e-6 * gamma/(2*pi) * B0; % kHz

% EPG simulations
npulse = 200;
alpha = alpha;
phi = mod(RF_phase_cycle(npulse, phi0), 2*pi);

repetitions = 1;

% single pool
[s0, Fn0, Zn0] = EPG_GRE(alpha*ones(npulse, 1), phi, TR, T1, T2);

% MT
tau = alpha/(gamma*B1);
energy = B1^2 * tau;
[smt, Fnmt, Znmt] = EPGX_GRE_MT(...
    alpha*ones(npulse, 1), phi, energy*ones(npulse, 1), ...
    TR, T1_MT, T2_MT, f_MT, k_MT, G);

% BM
[sx, Fnx, Znx] = EPGX_GRE_BM(...
    alpha*ones(npulse, 1), phi, TR, T1x, T2x, fx, kx, 'delta', delta_b);

save -v7 baseline/EPGX_GRE.mat Fn0 Zn0 Fnmt Znmt Fnx Znx
