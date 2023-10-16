clear
close all
%% Model parameters
m           = 0.2;      % kg
l           = 1;        % m
beta        = 0.003;    % kg m^2 / s
g           = 9.81;     % m / s^2

%% State space model 
A = [0 1;g/l -beta/(m*l^2)];
B = [0; 1/(m*l^2)];
%% LQR
Q = diag([1 1]);
R = 1;
[K, ~, ~] = lqr(A,B,Q,R);

%% Filter
% Using hint of (s+k)^2 or higher
k = 2;
lambda_1    = 2*k;              % 2*1*k
lambda_0    = (lambda_1/2)^2;   % k^2
% a = 5;
% b = 1;
% lambda_1 = a + b;
% lambda_0 = a*b;


%% Estimator paramters
% Gains 
g11         = 5000;             % 1 / ml^2
g12         = 0;                % cross gain?
g13         = 0;                % cross gain?
g21         = 0;                % cross gain? 
g22         = 5;                % beta / ml^2
g23         = 0;                % cross gain?
g31         = 0;                % cross gain? 
g32         = 0;                % cross gain? 
g33         = 50;               % 1 / l
Gamma       = [ g11 g12 g13;
                g21 g22 g23;
                g31 g32 g33];

m_0 = .1;
l_0 = 1.2;
beta_0 = 0.006;
% Initial estimates
theta_0 = [1/(m_0*l_0^2) beta_0/(m_0*l_0^2) 1/l_0]';           % 1 / ml^2 , beta / ml^2, 1 / l

% Normalization
alpha = 1;

% Input signals
A1 = 1;
A2 = .1;
w1 = pi;
w2 = pi*.3;

% Run Simulink
% dt              = 0.1;      % Sample time
t               = 15;        % Simulation time
sim("Problem_1.slx");

%% Plot
figure
sgtitle('Parameter estimates and their real values')
subplot(3,2,1)
hold on
plot(ans.actual.Time, ans.actual.Data(:,1))
plot([ans.actual.Time(1), ans.actual.Time(end)],[m m],'--')
hold off
ylabel('m')
legend('estimate','true value')
subplot(3,2,2)
hold on
plot(ans.theta.Time, ans.theta.Data(:,1))
plot([ans.theta.Time(1), ans.theta.Time(end)],[1/(m*l^2) 1/(m*l^2)],'--')
hold off
ylabel('1/(ml^2)')
legend('estimate','true value')
subplot(3,2,3)
hold on
plot(ans.actual.Time, ans.actual.Data(:,2))
plot([ans.actual.Time(1), ans.actual.Time(end)],[l l],'--')
hold off
ylabel('l')
legend('estimate','true value')
subplot(3,2,4)
hold on
plot(ans.theta.Time, ans.theta.Data(:,3))
plot([ans.theta.Time(1), ans.theta.Time(end)],[1/l 1/l],'--')
hold off
ylabel('l^{-1}')
legend('estimate','true value')
subplot(3,2,5)
hold on
plot(ans.actual.Time, ans.actual.Data(:,3))
plot([ans.actual.Time(1), ans.actual.Time(end)],[beta beta],'--')
hold off
ylabel('\beta')
legend('estimate','true value')
xlabel('t [s]')
subplot(3,2,6)
hold on
plot(ans.theta.Time, ans.theta.Data(:,2))
plot([ans.theta.Time(1), ans.theta.Time(end)],[beta/(m*l^2) beta/(m*l^2)],'--')
hold off
ylabel('\beta/(ml^2)')
legend('estimate','true value')
xlabel('t [s]')