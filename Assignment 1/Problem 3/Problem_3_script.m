%% SIM param
dt = 0.01;      % Sample time
%% Parameters
k_c = 2;        %Steady state gain
w_s = 10;       %Natural frequenzy
xi_s = 0.1;     %Damping ratio

w_h = .1;        % High Pass cutoff frequenzy
w_l = 1;        % Low Pass cutoff frequenzy
K = 10;         % Gain

a = .1;         % Sine wave amplitude
w = 3;          % Sine wave frequenzy

theta_0 = 0;    % Initial condition
%% Run simulink
sim("Problem_3_Simulink.slx");