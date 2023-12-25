clear variables;
clc;

% Simulation of the original conditions obtained by the LQR controller to
% the original non-linear system

% The system will have single derivative of the of the values and all state
% variables contribute to y_0

% Given values of constraints of the system
M = 1000; % Mass of the Crane in Kgs
m1 = 100; % mass of load 1 in Kgs
m2 = 100; % mass of load 2 in Kgs
l1 = 20; % length of cable 1 in meters
l2 = 10; % length of cable 2 in meters
g = 9.81;

% Initial conditions
y_0 = [5; 0; 5; 0; 10; 0];

% Defining the duration of the simulation
tspan = 0:0.01:9950;

% System parameters
params.M = M;
params.m1 = m1;
params.m2 = m2;
params.l1 = l1;
params.l2 = l2;
params.g = g;

%using inbuilt MATLAB function(ode45) to define the diff eqn
[t, y] = ode45(@(t, y) twoload(t, y, params), tspan, y_0);

% Plotting the function output on a 2D graph
figure
plot(t, y)
grid on
title('System Response')
xlabel('Time')
ylabel('State Variables')



function dydt = twoload(~, y, params)
    M = params.M;
    m1 = params.m1;
    m2 = params.m2;
    l1 = params.l1;
    l2 = params.l2;
    g = params.g;

    % Define system matrices
    A = [0 1 0 0 0 0;
        0 0 -m1*g/M 0 -m2*g/M 0;
        0 0 0 1 0 0;
        0 0 -((M+m1)*g)/(M*l1) 0 -m2*g/(M*l1) 0;
        0 0 0 0 0 1;
        0 0 -m1*g/(M*l2) 0 -(M*g+m2*g)/(M*l2) 0];

    B = [0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];
    
    % Defining Q and R matrices for LQR
    Q = [1500 0 0 0 0 0;
     0 1500 0 0 0 0;
     0 0 1500 0 0 0;
     0 0 0 1500 0 0;
     0 0 0 0 1500 0;
     0 0 0 0 0 1500];
    R=0.01;

   % Computing the LQR gain matrix
   K_Gainmat = lqr(A,B,Q,R);
   
   % control force based on the gain matrix and the current state
   F=-K_Gainmat*y;

   
    % Define system dynamics
    dydt = zeros(6, 1);
    dydt(1) = y(2);
    dydt(2)=(F-(g/2)*(m1*sind(2*y(3))+m2*sind(2*y(5)))-(m1*l1*(y(4)^2)*sind(y(3)))-(m2*l2*(y(6)^2)*sind(y(5))))/(M+m1*((sind(y(3)))^2)+m2*((sind(y(5)))^2));
    dydt(3) = y(4);
    dydt(4) = (dydt(2) * cosd(y(3)) - g * sind(y(3))) / l1;
    dydt(5) = y(6);
    dydt(6) = (dydt(2) * cosd(y(5)) - g * sind(y(5))) / l2;
end

