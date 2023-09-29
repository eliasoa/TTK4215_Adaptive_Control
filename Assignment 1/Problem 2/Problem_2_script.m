%% SIM param
dt = 0.01;  % Sample time

%% Parameters
theta_0 = 0;     % Initial condition

w_h = 1;    % High Pass cutoff frequenzy
w_l = 10;   % Low Pass cutoff frequenzy
K = 25;     % Gain

a = 0.1;    % Sine wave amplitude
w = 5;      % Sine wave frequenzy

%% Run simulink
sim("Problem_2_Simulink.slx");