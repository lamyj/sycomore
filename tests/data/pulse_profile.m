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

sampling_support_size = 501;

pulse_support = zero_crossings .* linspace(-1, 1, pulse_support_size + 1);
small_flip_angles = sinc(pulse_support);
small_flip_angles = small_flip_angles .* flip_angle / sum(small_flip_angles);
bandwidth = 2 * zero_crossings / pulse_duration;
slice_selection_gradient_moment = 2*pi*bandwidth*pulse_duration/slice_thickness;

sampling_locations = zeros(3, sampling_support_size);
sampling_locations(3, :) = linspace(...
    -slice_thickness, slice_thickness, sampling_support_size);

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

before_refocalization = zeros(3, size(sampling_locations, 2));
after_refocalization = zeros(3, size(sampling_locations, 2));

model.init_configuration(m0);

% Pulse
phase = pi;
model.RF(small_flip_angles(1), phase);
for j = 2 : numel(small_flip_angles)
    model.time(intervals.rf.mu, 'tau', intervals.rf.tau, 'p', intervals.rf.p);
    model.RF(small_flip_angles(j), phase);
end

% Slice profile before refocalization
for i=1:size(sampling_locations, 2)
    signal = model.isochromat(0, sampling_locations(:, i), []);
    before_refocalization(1, i) = real(signal.xy);
    before_refocalization(2, i) = imag(signal.xy);
    before_refocalization(3, i) = real(signal.z);
end

model.time(...
    intervals.refocalization.mu, 'tau', intervals.refocalization.tau, ...
    'p', intervals.refocalization.p);

% Slice profile after refocalization
for i=1:size(sampling_locations, 2)
    signal = model.isochromat(0, sampling_locations(:, i), []);
    after_refocalization(1, i) = real(signal.xy);
    after_refocalization(2, i) = imag(signal.xy);
    after_refocalization(3, i) = real(signal.z);
end

subplot(1, 3, 1);
plot(...
     pulse_duration .* linspace(-0.5, 0.5, pulse_support_size+1), ...
     small_flip_angles ./ max(small_flip_angles));
xlim([-0.5*pulse_duration, 0.5*pulse_duration]);
title('B^+_1(t)');

subplot(1, 3, 2);
plot(...
    sampling_locations(3, :), before_refocalization(1, :), ...
    sampling_locations(3, :), before_refocalization(2, :), ...
    sampling_locations(3, :), before_refocalization(3, :));
legend('m_x', 'm_y', 'm_z');
title('Before refocalization');

subplot(1, 3, 3);
plot(...
    sampling_locations(3, :), after_refocalization(1, :), ...
    sampling_locations(3, :), after_refocalization(2, :), ...
    sampling_locations(3, :), after_refocalization(3, :));
legend('m_x', 'm_y', 'm_z');
title('After refocalization');

fd = fopen('./baseline/pulse_profile.dat', 'w');
fwrite(fd, before_refocalization, 'float64');
fwrite(fd, after_refocalization, 'float64');
fclose(fd);
