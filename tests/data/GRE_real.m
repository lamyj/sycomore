T1 = 1000; % [ms]
T2 = 100; % [ms]
D = 0.89; % [µm^2/ms]
m0 = [0; 0; 1];

flip_angle = deg2rad(40); % [rad]
pulse_duration = 1; %[ms]
TR = 500; % [ms]
slice_thickness = 1000; % [um]
pulse_support_size = 100;
zero_crossings = 2;

TR_count = 10;

pulse_support = zero_crossings .* linspace(-1, 1, pulse_support_size + 1);
small_flip_angles = sinc(pulse_support);
small_flip_angles = small_flip_angles .* flip_angle / sum(small_flip_angles);
bandwidth = 2 * zero_crossings / pulse_duration;
slice_selection_gradient_moment = 2*pi*bandwidth*pulse_duration/slice_thickness;

model = CoMoTk;
model.R1 = 1/T1;
model.R2 = 1/T2;
model.D = D;

intervals = struct();

intervals.rf.mu = 1+numel(fieldnames(intervals));
intervals.rf.tau = pulse_duration/pulse_support_size;
intervals.rf.p = [0; 0; slice_selection_gradient_moment / pulse_support_size];

intervals.echo.mu = 1+numel(fieldnames(intervals));
intervals.echo.tau = (TR - pulse_duration)/2.;
intervals.echo.p = [0; 0; -slice_selection_gradient_moment / 2.];

signal = zeros(3, TR_count);
model.init_configuration(m0);

% % time after echo
% cm_bSSFP_real.time( mu_bSSFP_real, 'tau', tau_bSSFP_real, 'p', p_bSSFP_real );
% % excitation pulse
% cm_bSSFP_real.RF( al_rad( 1 ), ph_rad( i ) );     
% for j = 1 : n_tau
%     cm_bSSFP_real.time( mu_real_rf, 'tau', tau_real_rf, 'p', p_real_rf );
%     cm_bSSFP_real.RF( al_rad( j + 1 ), ph_rad( i ) );     
% end
% cm_bSSFP_real.time( mu_bSSFP_real, 'tau', tau_bSSFP_real, 'p', p_bSSFP_real );
% % bSSFP : only the slice encoding direction requires a zero gradient moment
% b_n = cm_bSSFP_real.b_n & reshape( cm_bSSFP_real.p_n( 3, : ) == 0, size( cm_bSSFP_real.b_n ) );
% iso = cm_bSSFP_real.isochromat( 0, [], b_n );

tic;
for i = 1 : TR_count
    % Pulse
    phase = pi/3+mod(i-1, 2)*pi;
    model.RF(small_flip_angles(1), phase); 
    for j = 2 : numel(small_flip_angles)
        model.time(intervals.rf.mu, 'tau', intervals.rf.tau, 'p', intervals.rf.p);
        model.RF(small_flip_angles(j), phase);     
    end
    
    model.time(intervals.echo.mu, 'tau', intervals.echo.tau, 'p', intervals.echo.p);
    
    % WARNING: which configurations contribute to the signal intensity?
    s = model.isochromat(0, [], []);
    signal(:,i) = [real(s.xy), imag(s.xy), real(s.z)];
    
    model.time(intervals.echo.mu, 'tau', intervals.echo.tau, 'p', intervals.echo.p);
end
toc;

fd = fopen('./baseline/GRE_real.dat', 'w');
fwrite(fd, signal, 'float64');
fclose(fd);