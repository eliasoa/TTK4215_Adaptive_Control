close 
clear
%% Plant
b0                  = .5;
b1                  = .3;
a0                  = .3;
a1                  = .2;
kp                  = b1;
Zp                  = [0 1 b0/b1];
Rp                  = [1 a1 a0];
zp                  = tf(kp*Zp,Rp,-1)

%% Model
km                  = 1;
Zm                  = 1;
Rm                  = [1 0.1];
%% Filters
lambda0             = .1;
Lambda0             = [1 lambda0];
L_z                 = tf(1,Lambda0)
%% Actual params
S                   = [ 1  0   b1  0;
                        a1  b1  b1*lambda0 + b0 0;
                        a0  b0  lambda0*b0  0;
                        0   0   0   1];
p                   = [ a1 - 0.1 - b0/b1;
                        a0 + a1*lambda0 - 0.1*(b0/b1) - 0.1*lambda0 - (b0/b1)*lambda0;
                        lambda0*a0 - 0.1*(b0/b1)*lambda0;
                        km/kp];
theta_bar           = S\p;
%% Estimator
gamma               = 1;
theta_0             = [1 1 1 1]';    % Initial guess for controller params
alpha               = 1;
%% Simulation
r1                  = 1.2;
freq1               = .2;
r2                  = 1;
freq2               = .23;
r3                  = .5;
freq3               = .5;
simTime             = 500;   
sim('BRGH_MRAC_simulink.slx')