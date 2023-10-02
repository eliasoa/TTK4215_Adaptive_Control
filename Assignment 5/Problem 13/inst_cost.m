close all;
clear;
%% Insert true system
beta                    = 0.2;
m                       = 15;
k                       = 2;
A                       = [0 0 0; 0 0 1;0 0 -beta/m];               
B                       = [0; 0; 1/m];
%% Define filters
lambda_1                = .1; 
lambda_0                = 6;
DEN                     = [1 lambda_1 lambda_0];
[  ~,  ~,C_f_1,D_f_1]   = tf2ss([0 0 1],DEN);     
[  ~,  ~,C_f_2,D_f_2]   = tf2ss([0 1 0],DEN);                            
[A_f,B_f,C_f_3,D_f_3]   = tf2ss([1 0 0],DEN);                     
% A_f and B_f is only needed once as it is the same all over. Hence tildes

%% Simulation MATLAB
h           = 0.01;  % sample time (s)
N           = 1200; % number of samples

t           = 0:h:h*(N-1);

% Define input as a function of t
u           = 5 * sin(2 * t) + 10.5;
gamma_1     = 1;        
gamma_2     = diag([50 1]);

% Memory allocation
x           = zeros(3, N);
x_z_1       = zeros(1,1);
x_z_2       = zeros(2, 1);
x_phi_1     = zeros(1, 1);
x_phi_2     = zeros(2, 1);
theta_1     = zeros(1, N);  % estimate of k
theta_2     = zeros(2, N);  % estimate of m and beta


% Initial estimates
theta_1(:,1) = 1;       % Initial estimate of k
theta_2(:,1) = [1; 1];   % Initial estimate of m and beta
%% Main loop. Simulate using forward Euler (x[k+1] = x[k] + h*x_dot)
for n = 1:N-1
 
    % Simulate true system
    x_dot           = A*x(:, n) + B*u(n);
    x(:, n+1)       = x(:, n) + h*x_dot;
    y_2             = x(2, n);
    y_1             = u(n)/k + y_2;

    % Generate z and phi by filtering known signals
    x_z_2_n           = x_z_2 + (A_f*x_z_2 + B_f*u(n))*h; 
    z_2               = C_f_1*x_z_2;  % z_2 = 1/Lambda [u]
                    
    
    x_phi_2_n         = x_phi_2 + (A_f*x_phi_2 + B_f*y_2)*h;
                                            
    phi_2             = [ (C_f_3*x_phi_2 + D_f_3*y_2);  % s^2/Lambda * y_2
                          (C_f_2*x_phi_2 + D_f_2*y_2)   % s/Lamnba * y2
                        ];
    
    % Generate first system
    z_1 = u(n);
    phi_1 = y_1 - y_2;
    
    % Calculate estimation error
    n_s_1 = phi_1'*phi_1;
    n_s_2 = phi_2'*phi_2;

    epsilon_1           = z_1 - theta_1(:,n)' * phi_1; 
    epsilon_2           = z_2 - theta_2(:,n)' * phi_2;

    % Update law
    theta_2_dot         = gamma_2 * epsilon_2 * phi_2;
    theta_2(:, n+1)     = theta_2(:, n) + theta_2_dot * h;
    
    theta_1_dot = gamma_1*epsilon_1*phi_1;
    theta_1(:, n+1) = theta_1(:, n) + theta_1_dot*h;

    % Set values for next iteration
    x_phi_2           = x_phi_2_n;
    x_z_2             = x_z_2_n;
end

% Plots
figure
subplot(3,1,1)
plot(t, theta_2(1,:)); hold on
plot([t(1), t(end)],[m, m]);
hold off
ylabel('m')
legend('estimate','true value')
grid
subplot(3,1,2)
plot(t, theta_2(2,:)); hold on
plot([t(1), t(end)],[beta, beta]); 
hold off
ylabel('\beta')
legend('estimate','true value')
grid
subplot(3,1,3)
plot(t, theta_1(1,:)); hold on
plot([t(1), t(end)],[k, k]); 
hold off
ylabel('k')
grid
legend('estimate','true value')
title('Parameter estimates')
xlabel('t [s]')