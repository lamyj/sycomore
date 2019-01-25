T1 = inf; % [ms]
T2 = inf; % [ms]
D = 0; % [um^2/ms]
m0 = [0; 0; 1];

flip_angle = deg2rad(90); % [rad]
pulse_duration = 1; % [ms]
% NOTE: in the absence of relaxation and diffusion, the TR is meaningless
TR = 500; % [ms]
slice_thickness = 1000; % [um]
pulse_support_size = 100;
zero_crossings = 2;

pulse_support = zero_crossings .* linspace(-1, 1, pulse_support_size + 1);
small_flip_angles = sinc(pulse_support);
small_flip_angles = small_flip_angles .* flip_angle / sum(small_flip_angles);
bandwidth = 2 * zero_crossings / pulse_duration;
slice_selection_gradient_moment = 2*pi*bandwidth*pulse_duration/slice_thickness;

frequencies = linspace(-30, 30, 201);

model = CoMoTk;
model.R1 = 1/T1;
model.R2 = 1/T2;
model.D = D;

intervals = struct();

intervals.rf.mu = 1+numel(fieldnames(intervals));
intervals.rf.tau = pulse_duration/pulse_support_size;
intervals.rf.p = [0; 0; slice_selection_gradient_moment / pulse_support_size];

intervals.refocalization.mu = 1+numel(fieldnames(intervals));
intervals.refocalization.tau = (TR - pulse_duration)/2.;
intervals.refocalization.p = [0; 0; -slice_selection_gradient_moment / 2];

model.init_configuration(m0);

% Pulse
phase = pi;
model.RF(small_flip_angles(1), phase);
for j = 2 : numel(small_flip_angles)
    model.time(intervals.rf.mu, 'tau', intervals.rf.tau, 'p', intervals.rf.p);
    model.RF(small_flip_angles(j), phase);
end

model.time(...
    intervals.refocalization.mu, 'tau', intervals.refocalization.tau, ...
    'p', intervals.refocalization.p);

% Slice profile after refocalization
magnetization = zeros(2, numel(frequencies));
for i=1:numel(frequencies)
    signal = model.isochromat(frequencies(i), [], []);
    magnetization(1, i) = abs(signal.xy);
    magnetization(2, i) = real(signal.z);
end

subplot(1, 1, 1);
plot(...
    frequencies, magnetization(1,:), ...
    frequencies, magnetization(2,:));
legend('Transversal', 'Longitudinal');
xlabel('Relative frequency (rad/ms)');
title('Magnetization');

fd = fopen('./baseline/off_resonance.dat', 'w');
fwrite(fd, magnetization, 'float64');
fclose(fd);
