clear; clc;

% FIRST ORDER WITH DECAY --------------------------------------------------
% Stimulus = mRNA availability (Binary)
% Response = Protein synthesis

alpha = 0.8;  % Decay constant
beta = 1.0;   % mRNA switch
t1 = 2; 
t2 = 8; 
tSpan = [0, 15];
z0 = 0;

% Stimulus (mRNA formation)
stimulus = @(t) beta * (t > t1 & t < t2);

% System of ODEs
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
dzdt = @(t, z) stimulus(t) - alpha * z;
[t, z] = ode45(dzdt, tSpan, z0, options);

% Normalize the response
z_norm = z / max(z);

stimulus_arr = arrayfun(stimulus, t);

figure;
ylim([0 1.2]);
hold on;
plot(t, stimulus_arr, 'r', 'DisplayName', 'Stimulus (mRNA)');
plot(t, z_norm, 'b', 'DisplayName', 'Response (Protein conc.)');
xlabel('Time');
ylabel('Normalized signal & response');
title('First order dynamics with decay');
legend('Location', 'best');
grid on;
hold off;
