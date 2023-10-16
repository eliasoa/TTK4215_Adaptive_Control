close all
clear
%% Model parameters
% Reference model
am      = 2;
bm      = 2;

% Real model
a       = 1;
b       = 5;        % Unknown

%% Parameter estimator
gamma1  = 10;     %k
gamma2  = 10;    %l

%% Other
ref     = 1;
Ts      = 5;        % Time constant reference signal
a       = 1;        % Amplitude time varying reference signal
w       = 1;        % Frequency --||--

%% Run simulation
t       = 50;       % Simulation time

sim("Problem2.slx");