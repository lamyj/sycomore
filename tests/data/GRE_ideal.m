T1 = 1000; % [ms]
T2 = 100; % [ms]
D = 0.89; % [um^2/ms]
m0 = [0; 0; 1];

flip_angle = deg2rad(40); % [rad]
TR = 500; % [ms]

TR_count = 10;

model = CoMoTk;
model.R1 = 1/T1;
model.R2 = 1/T2;
model.D = D;

intervals = struct();

intervals.echo.mu = 1+numel(fieldnames(intervals));
intervals.echo.tau = TR/2;

signal = zeros(3, TR_count);
model.init_configuration(m0);

for i = 1 : TR_count
    model.RF(flip_angle, pi/3+mod(i-1, 2)*pi);
    model.time(intervals.echo.mu, 'tau', intervals.echo.tau);
    s = model.isochromat(0, [], []);
    signal(:,i) = [real(s.xy), imag(s.xy), real(s.z)];
    model.time(intervals.echo.mu, 'tau', intervals.echo.tau);
end

fd = fopen('./baseline/GRE_ideal.dat', 'w');
fwrite(fd, signal, 'float64');
fclose(fd);