clc; clear;

% Complex Network of Feed Forward Loops -----------------------------------
% X1 is activated by environmental stimulus (Binary)
% X1 --> activates --> Y1
% X1, Y1 --> Incoherent FFL --> Z1
% X1, Y1 --> Coherent FFL --> X2
% X2 --> activates --> Y2
% X2, Y2 --> Incoherent FFL --> Z2
% X2, Y2 --> Coherent FFL --> Z3

alphay1 = 0.4;
betay1 = 1.5;

alphaz1 = 0.45;
betaz1 = 1.3;

alphax2 = 0.5;
betax2 = 1.25;

alphay2 = 0.55;
betay2 = 1.15;

alphaz2 = 0.6;
betaz2 = 1.2;

alphaz3 = 0.65;
betaz3 = 1.0;

tSpan = [0, 20];

Kx1z1 = 0.5;
Ky1z1 = 0.6;
Kx2y2 = 0.4;

Kx2z2 = 0.1;
Ky2z2 = 0.8;

Kx2z3 = 0.2132;
Ky2z3 = 0.4855;

% initial conditions
y10 = 0;
z10 = 0;
x20 = 0;
y20 = 0;
z20 = 0;
z30 = 0;

activationTime = 2;
deactivationTime = 8;

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

t = 0:0.01:max(tSpan);
x1 = @(t) 1 * (t > activationTime & t < deactivationTime);
x_arr = arrayfun(x1, t);

dy1_dt = @(t1, y1) (betay1 * (x1(t1) == 1)) - alphay1 * y1;
[t1, y1] = ode45(dy1_dt, tSpan, y10, options);

y1_interpolated = @(t2) interp1(t1, y1, t2, 'linear', 'extrap');

dz1_dt = @(t2, z1) (betaz1 * (x1(t2) == 1 & y1_interpolated(t2) < Kx1z1 * max(y1)) - alphaz1 * z1);
[t2, z1] = ode45(dz1_dt, tSpan, z10, options);

dx2_dt = @(t3, x2) (betax2 * (x1(t3) == 1 & y1_interpolated(t3) > Ky1z1) - alphax2 * x2);
[t3, x2] = ode45(dx2_dt, tSpan, x20, options);

x2_interpolated = @(t4) interp1(t3, x2, t4, 'linear', 'extrap');

dy2_dt = @(t4, y2) (betay2 * (x2_interpolated(t4) > Kx2y2 * max(y2))) - alphay2 * y2;
[t4, y2] = ode45(dy2_dt, tSpan, y20, options);

y2_interpolated = @(t5) interp1(t4, y2, t5, 'linear', 'extrap');
x2_interpolated = @(t5) interp1(t3, x2, t5, 'linear', 'extrap');

dz2_dt = @(t5, z2) ((betaz2 * ((x2_interpolated(t5) > Kx2z2 * max(x2)) & (y2_interpolated(t5) < Ky2z2 * max(y2)))) - alphaz2 * z2);
[t5, z2] = ode45(dz2_dt, tSpan, z20, options);

dz3_dt = @(t6, z3) ((betaz3 * ((x2_interpolated(t6) > Kx2z3 * max(x2)) & (y2_interpolated(t6) > Ky2z3 * max(y2)))) - alphaz3 * z3);
[t6, z3] = ode45(dz3_dt, tSpan, z30, options);

figure;
hold on;
plot(t, x_arr, 'DisplayName', 'X1 [Binary]');
plot(t1, y1, 'DisplayName', 'Y1');
plot(t2, z1, 'DisplayName', 'Z1', 'LineWidth', 1.1);
plot(t3, x2, 'DisplayName', 'X2', 'LineWidth', 1.2);
plot(t4, y2, 'DisplayName', 'Y2', 'LineWidth', 1.3);
plot(t5, z2, 'DisplayName', 'Z2', 'LineWidth', 1.4);
plot(t6, z3, 'DisplayName', 'Z3', 'LineWidth', 1.5);

xlabel('Time');
ylabel('Gene transcription dynamics');
title('Complex FFL circuit');
legend('Location', 'northeast');
grid on;
hold off;