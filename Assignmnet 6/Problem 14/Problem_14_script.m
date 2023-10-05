%% Model parameters
[k0 w0 xi0 k1 w1 xi1] = linearInterpolation(20);

theta_1_real = [k0*w0^2 2*xi0*w0 w0^2]
theta_2_real = [k1*w1^2 2*xi1*w1 w1^2]

%% Filter
lambda_1 = 20;
lambda_0 = 100;

%% Estimator paramters
% Gains 
Gamma_1 = diag([7 8 10]);    %  FORTEGNSFEIL ET STED??? 
Gamma_2 = diag([1 1 1]);
theta_1_0 = [200 9 350]';       % kw^2 2*xi*w w^2
theta_2_0 = [10 8 200]';       % kw^2 2*xi*w w^2

% No normalization
alpha = 0;

%% Run Simulink
dt              = 0.01;      % Sample time
t               = 500;        % Simulation time
sim("Problem_14.slx")

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