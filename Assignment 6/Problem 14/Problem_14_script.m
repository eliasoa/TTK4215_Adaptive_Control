% clear all
% close all
%% Model parameters
[k0 w0 xi0 k1 w1 xi1] = linearInterpolation(20);

theta_1_real = [k0*w0^2 2*xi0*w0 w0^2]
theta_2_real = [k1*w1^2 2*xi1*w1 w1^2]

%% Filter
lambda_1 = 20;
lambda_0 = 100;

%% Estimator paramters
% Gains 
Gamma_1 = diag([1000 -5 -100]);       % kw^2 2*xi*w w^2
Gamma_2 = diag([20 100000 1200]);

% Use know values at V = 30 as initial conditions
[k0_init w0_init xi0_init k1_init w1_init xi1_init] = linearInterpolation(30);
theta_1_0 = [k0_init*w0_init^2 2*xi0_init*w0_init w0_init^2]'% kw^2 2*xi*w w^2
theta_2_0 = [k1_init*w1_init^2 2*xi1_init*w1_init w1_init^2]'% kw^2 2*xi*w w^2

% Normalization
alpha = 1;
% Run Simulink
dt              = 0.01;      % Sample time
t               = 150;        % Simulation time
sim("Problem_14.slx");

%% Plot
figure
sgtitle('Parameter estimates and their real values')
subplot(3,2,1)
hold on
plot(ans.theta_1.Time, ans.theta_1.Data(:,1))
plot([ans.theta_1.Time(1), ans.theta_1.Time(end)],[w0 w0],'--')
hold off
ylabel('\omega_0')
legend('estimate','true value')
subplot(3,2,2)
hold on
plot(ans.theta_2.Time, ans.theta_2.Data(:,1))
plot([ans.theta_2.Time(1), ans.theta_2.Time(end)],[w1 w1],'--')
hold off
ylabel('\omega_1')
legend('estimate','true value')
subplot(3,2,3)
hold on
plot(ans.theta_1.Time, ans.theta_1.Data(:,2))
plot([ans.theta_1.Time(1), ans.theta_1.Time(end)],[xi0 xi0],'--')
hold off
ylabel('\xi_0')
legend('estimate','true value')
subplot(3,2,4)
hold on
plot(ans.theta_2.Time, ans.theta_2.Data(:,2))
plot([ans.theta_2.Time(1), ans.theta_2.Time(end)],[xi1 xi1],'--')
hold off
ylabel('\xi_1')
legend('estimate','true value')
subplot(3,2,5)
hold on
plot(ans.theta_1.Time, ans.theta_1.Data(:,3))
plot([ans.theta_1.Time(1), ans.theta_1.Time(end)],[k0 k0],'--')
hold off
ylabel('k_0')
legend('estimate','true value')
xlabel('t [s]')
subplot(3,2,6)
hold on
plot(ans.theta_2.Time, ans.theta_2.Data(:,3))
plot([ans.theta_2.Time(1), ans.theta_2.Time(end)],[k1 k1],'--')
hold off
ylabel('k_1')
legend('estimate','true value')
xlabel('t [s]')
%% Functions
function [k0 w0 xi0 k1 w1 xi1] = linearInterpolation(V)
    % v    k0   w0   xi0  k1    w1   xi1
v=  [30  0.81 19.75 0.31 0.064 14.0 0.365; 
     60  0.77 19.0  0.27 0.09  13.5 0.505];
a = zeros(6,1);
for i = 1:size(a)
    a(i) = (v(2,i+1) - v(1,i+1))/(v(2,1)-v(1,1));
end
b = zeros(6,1);
for i = 1:size(b)
    b(i) = v(1,i+1) - a(i)*v(1,1);
end

param = zeros(6);
for i = 1:6
    param(i) = a(i)*V +b(i);
end
k0  = param(1);
w0  = param(2);
xi0 = param(3);
k1  = param(4);
w1  = param(5);
xi1 = param(6);
end