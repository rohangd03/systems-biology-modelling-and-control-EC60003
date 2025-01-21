clc; clear;

% Declare global variables
global t1 t2 t3 E V1 V2 V3 I

% Parameter initialization
t1 = 6;
t2 = 100;
t3 = 36;
E = 0.2;
V1 = 3;
V2 = 11;
V3 = 10;
I = 216;

% Simulation settings
simulationTime = 1000;
dt = 0.5;
t = 0:dt:simulationTime; % Time vector

steps = length(t); % Total number of steps

% Preallocate variables
x = zeros(steps, 1);
y = zeros(steps, 1);
z = zeros(steps, 1);
h1 = zeros(steps, 1);
h2 = zeros(steps, 1);
h3 = zeros(steps, 1);

for i = 1:steps-1
    global t1 t2 t3 E V1 V2 V3 I

    [f1, f2, f3, f4, f5] = calc(z(i), y(i), h3(i));

    
    dxdt = f1 - E * (x(i) / V1 - y(i) / V2) - x(i) / t1;
    dydt = E * (x(i) / V1 - y(i) / V2) - y(i) / t2;
    dzdt = f5 + I - f2 - f3 * f4;
    dh1dt = 3 * (x(i) - h1(i)) / t3;
    dh2dt = 3 * (h1(i) - h2(i)) / t3;
    dh3dt = 3 * (h2(i) - h3(i)) / t3;

    
    x(i+1) = x(i) + dxdt * dt;
    y(i+1) = y(i) + dydt * dt;
    z(i+1) = z(i) + dzdt * dt;
    h1(i+1) = h1(i) + dh1dt * dt;
    h2(i+1) = h2(i) + dh2dt * dt;
    h3(i+1) = h3(i) + dh3dt * dt;
end

figure;
plot(t, x);
hold on;
title('x vs Time');
xlabel('Time (min)');
ylabel('x');
grid on;

figure;
plot(t, y);
title('y vs Time');
xlabel('Time (min)');
ylabel('y');
grid on;

figure;
plot(t, z);
title('Glucose vs Time');
xlabel('Time (min)');
ylabel('z');
grid on;

figure;
plot(t, h1);
title('h1 vs Time');
xlabel('Time (min)');
ylabel('h1');
grid on;

figure;
plot(t, h2);
title('h2 vs Time');
xlabel('Time (min)');
ylabel('h2');
grid on;

figure;
plot(t, h3);
title('h3 vs Time');
xlabel('Time (min)');
ylabel('h3');
grid on;

function [f1, f2, f3, f4, f5] = calc(z, y, h3)
    global t2 E V1 V2 V3

    % Compute the intermediate functions
    f1 = 209 / (1 + exp(-z / (300 * V3) + 6.6));
    f2 = 72 * (1 - exp(-z / (144 * V3)));
    f3 = 0.01 * z / V3;
    f4 = 90 / (1 + exp(-1.772 * log((y / V2) + (1 / (E*t2))) + 7.76)) + 4;
    f5 = 180 / (1 + exp(((0.29 * h3) / V1) - 7.5));
end
