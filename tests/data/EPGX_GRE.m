% EPG-X, revision 5da0e7128d2c7fc42ff13991ddeba806c81ed23f

TR = 5; % ms
alpha = deg2rad(10); % rad
phi0 = 117; % deg

% Single pool
T1=779; % ms
T2=45; % ms

% MT
T1_MT = [779 779]; % ms
f_MT = 0.117; % unitless
k_MT = 4.3e-3; % kHz
T2_MT = 45; % ms

% Exchange
T1x = [1000 500]; % ms
T2x = [100 20]; % ms
kx = 2e-3; % kHz
fx = 0.2; % unitless

% RF saturation factor for MT
G = 15.1; % us 
B1 = 13; % uT
gamma = 267.5221 * 1e-3; % rad/ms/uT
t_RF = alpha/(gamma*B1); % ms
b1sqrdtau = B1^2 * t_RF;

% EPG simulations
npulse = 200;
alpha = alpha*ones(npulse, 1);
phi = mod(RF_phase_cycle(npulse, phi0), 2*pi);

repetitions = 1;

% single pool
[s0, Fn0, Zn0] = EPG_GRE(alpha, phi, TR, T1, T2);

% MT
[smt, Fnmt, Znmt] = EPGX_GRE_MT(...
    alpha, phi, b1sqrdtau*ones(npulse, 1), ...
    TR, T1_MT, T2_MT, f_MT, k_MT, G);

% BM
[sx, Fnx, Znx] = EPGX_GRE_BM(alpha, phi, TR, T1x, T2x, fx, kx);

save -v7 baseline/EPGX_GRE.mat Fn0 Zn0 Fnmt Znmt Fnx Znx
