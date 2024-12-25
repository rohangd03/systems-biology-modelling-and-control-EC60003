clear; clc;

% --------------------------------------------------------------------
% [1] SIR (Base model)

% Simulation Parameters
% beta = Infection rate
% gamma = Recovery rate

N = 10000; % Total number of people in population
infection_rate = (0.002/7)*N;  
recovery_rate = 0.15;
simulationTime = 50;

t1 = 0:0.1:simulationTime;

S_initial = 0.99;
I_initial = 0.01;
R_initial = 0;
initial = [S_initial, I_initial, R_initial];
% ---------------------------------------------------------------
[S, I, R] = sir_model(t1, infection_rate, recovery_rate, initial);

figure;
grid on;
plot(t1, S, 'linewidth', 1);
hold on;
plot(t1, I, 'linewidth', 1);
hold on;
plot(t1, R, 'linewidth', 1);
legend({'Susceptible', 'Infected', 'Recovered'});
ylabel('Fraction of population');
xlabel('time (days)');
title('SIR Model (Total population = 10000)');
grid on;
hold off;

% ----------------------------------------------------------------------
% [2] Variation of susceptibility with infection rate
% Reproduction number = Infection rate/Recovery rate

R_array = [0.002, 0.004, 0.006, 0.008, 0.01];
N = 10000;
S0 = 0.99 * N;           
I0 = 0.01 * N;             

figure;
colors = lines(length(R_array));
for i = 1:length(R_array)
    R = R_array(i);
    
    Susceptible = linspace(S0, 0.01 * N, 1000);

    dI_dS = -1 + 1 ./ (R * Susceptible);
    Infected = I0 + cumtrapz(Susceptible, dI_dS);
    
    normalizedS = Susceptible/N;
    normalizedI = Infected/N;
    
    plot(normalizedS, normalizedI, 'Color', colors(i, :), 'LineWidth', 1);
    hold on;
    legendStrings{i} = sprintf('R = %.3f', R);
end

xlabel('Normalized Susceptible population (S)');
ylabel('Normalized Infected population (I)');
title('Variability of I with S (Concave function)');
legend(legendStrings);
grid on;
hold off;

% ----------------------------------------------------------------------
% [3] SIR Model with LOSS OF IMMUNITY

% Simulation Parameters -----------------------------------------
% beta = Infection rate
% gamma = Recovery rate
% alpha = Rate of loss of immunity

N = 10000; % Total number of people in population
infection_rate = (0.002/7)*N;  
recovery_rate = 0.15;
immunity_loss = 0.05;

simulationTime = 50;

t = 0:0.1:simulationTime;

S_initial = 0.99;
I_initial = 0.01;
R_initial = 0;

initial = [S_initial, I_initial, R_initial];

[S, I, R] = sir_loss_of_immunity(t, infection_rate, recovery_rate, immunity_loss, initial);

figure;
grid on;
plot(t, S, 'linewidth', 1);
hold on;
plot(t, I, 'linewidth', 1);
hold on;
plot(t, R, 'linewidth', 1);
legend({'Susceptible', 'Infected', 'Recovered'});
ylabel('Fraction of population');
xlabel('time (days)');
title('SIR Model with loss of immunity');
grid on;
hold off;

% ---------------------------------------------------------------------
% [4] Endemic state vs Endemic free state is determined by Reproduction number

N = 10000; % Total number of people in population

t = 0:0.1:40;

S_initial = 0.99 * N;
I_initial = 0.01 * N;
R_initial = 0;

initial = [S_initial, I_initial, R_initial];
recovery_rate = 0.15; % fixed
immunity_loss = 0.005;

R_array = linspace(0, 2, 100);
I_array = zeros(length(R_array), 1);

for i = 1:length(R_array)
    R = R_array(i);
    infection_rate = ((recovery_rate + immunity_loss) * R)/N;
    [S, I, R] = sir_loss_of_immunity(t, infection_rate, recovery_rate, immunity_loss, initial);
    I_final = I(end);
    I_array(i) = I_final;
end

figure;
x = [1, 1];
ylim([0 I_array(end)+1000]);
y = ylim;
plot(R_array, I_array, 'r-', 'LineWidth', 1);
hold on;
plot(x, y, 'k:', 'LineWidth', 1.5); % Plot dotted line of x = 1
xlabel('Reproduction number');
ylabel('Steady state Infected population');
title('Bifurcation of I with R (Endemic free state --> Endemic state)');
grid on;

% ------------------------------------------------------------------------
% [5] SIR model with Vital dynamics (generic death)

% Simulation Parameters -----------------------------------------
% beta = Infection rate
% gamma = Recovery rate
% mu = Rate of generic death
% B = Total "dead" population
% S, I, R, B normalized w.r.t to total population
% S = y(1); I = y(2); R = y(3)

N = 10000; % Normalized population
infection_rate = (0.002/7) * N;  
recovery_rate = 0.15;
death_rate = (0.002/7) * 0.1;

simulationTime = 50;

t = 0:0.1:simulationTime;

S_initial = 0.99;
I_initial = 0.01;
R_initial = 0;

initial = [S_initial, I_initial, R_initial];
% ---------------------------------------------------------------
[S, I, R] = sir_vital_dynamics(t, infection_rate, recovery_rate, death_rate, initial);

figure;
grid on;
plot(t, S, 'linewidth', 1);
hold on;
plot(t, I, 'linewidth', 1);
hold on;
plot(t, R, 'linewidth', 1);

legend({'Susceptible', 'Infected', 'Recovered'});
ylabel('Fraction of population');
xlabel('time (days)');
title('SIR Model with Vital Dynamics');
grid on;
hold off;

% ------------------------------------------------------------------------
% --------------------- ALL THE FUNCTIONS --------------------------------

function [S,I,R] = sir_model(t, infection_rate, recovery_rate, initial)
    dydt = @(t,y) [
        (- infection_rate * y(1) * y(2)); 
        (infection_rate * y(1)*y(2) - recovery_rate * y(2)); 
        (recovery_rate * y(2))];

    [~, y] = ode45(dydt, t, initial);
    S = y(:,1);
    I = y(:,2);
    R = y(:,3);
end

function [S,I,R] = sir_loss_of_immunity(t, infection_rate, recovery_rate, immunity_loss, initial)
    dydt = @(t,y) [(- infection_rate * y(1) * y(2) + immunity_loss * y(3)); 
                   (infection_rate * y(1) * y(2) - recovery_rate * y(2)); 
                   (recovery_rate * y(2) - immunity_loss * y(3))];
    [~, y] = ode45(dydt, t, initial);
    S = y(:,1);
    I = y(:,2);
    R = y(:,3);
end

function [S,I,R] = sir_vital_dynamics(t, infection_rate, recovery_rate, death_rate, initial)
    dydt = @(t,y) [(- infection_rate * y(1) * y(2) + death_rate - death_rate * y(1));
                   (infection_rate * y(1) * y(2) - recovery_rate * y(2) - death_rate * y(2));
                   (recovery_rate * y(2) - death_rate * y(3))];
    [~,y] = ode45(dydt, t, initial);
    S = y(:,1);
    I = y(:,2);
    R = y(:,3);
end