clear
close all
%% Model parameters
m                       = 0.2;      % kg
l                       = 1;        % m
beta                    = 0.003;    % kg m^2 / s
g                       = 9.81;     % m / s^2

%% State space model 
A                       = [0 1;g/l -beta/(m*l^2)];
B                       = [0; 1/(m*l^2)];
%% Inital conditions and saturation
limit_q                 = 14.76; % limit at +- 14.76
limit_dot_q             = 0;
q_0                     = 0;    %deg2rad(0);  
q_dot_0                 = -0.01; %deg2rad(10);

sat_limit               = 0.5;  % N*m

%% LQR
Q                       = diag([1 1]);
R                       = 1;
[K, ~, ~]               = lqr(A,B,Q,R);

%% Lyapunov
k = 1;

%% Switching logic
angle_threshold = deg2rad(160);

%% Parameter Estimator
useEstimates            = 1;                % for probelm g
switch_time             = 15;               % ----||----

% Filter
k_f                     = 2;
lambda_1                = 2*k_f;              % 2*1*k
lambda_0                = (lambda_1/2)^2;       % k^2

% Initial estimates
m_0                     = .1;
l_0                     = 1.2;
beta_0                  = 0.006;
theta_0                 = [1/(m_0*l_0^2) beta_0/(m_0*l_0^2) 1/l_0]';           % 1 / ml^2 , beta / ml^2, 1 / l

% Adaptive gain
Gamma                   = diag([5000 5 50]);
% Normalization
alpha                   = 1;

% Input signals
A1                      = 0.95;
A2                      = .1;
w1                      = pi*.9;
w2                      = pi/10;

%% Run simulation
t                        = 10;          % Simulation time
sim("A8_Problem_1.slx");

%% Plot
simout = ans;
t_state = simout.states.Time;
q = simout.states.Data(:,1);
q_dot = simout.states.Data(:,2);

t_tau = simout.tau.Time;
unsat = simout.tau.Data(:,1);
sat = simout.tau.Data(:,2);

figure(1)
subplot(311)
hold on
plot(t_state, rad2deg(q));
plot([t_state(1) t_state(end)], [180 180],'--');
hold off
legend('q','ref')
ylabel('angle [deg]')

subplot(312)
hold on
plot(t_state, rad2deg(q_dot));
plot([t_state(1) t_state(end)], [0 0],'--');
hold off
legend('$\dot{q}$','ref','Interpreter','latex')
ylabel('angular velocity [deg/s]')

subplot(313)
hold on
plot(t_tau,unsat)
plot(t_tau,sat,'--')
hold off
xlabel('time [s]')
ylabel('\tau [Nm]')
legend('Commanded','Saturated')
