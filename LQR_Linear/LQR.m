% LQR Controller
clear variables;
clc;

% Given values of constraints of the system
M = 1000; % Mass of the Crane in Kgs
m1 = 100; % mass of load 1 in Kgs
m2 = 100; % mass of load 2 in Kgs
l1 = 20; % length of cable 1 in meters
l2 = 10; % length of cable 2 in meters
g = 9.81;

% Values of matrix A and B from the state-space representation of the linearized system
A = [0 1 0 0 0 0;
     0 0 -(g*m1)/M 0 -(g*m2)/M 0;
     0 0 0 1 0 0;
     0 0 (-g*M -m1*g)/(M*l1) 0 -(g*m2)/(M*l1) 0;
     0 0 0 0 0 1;
     0 0 -(g*m1)/(M*l2) 0 (-g*M -g*m2)/(M*l2) 0];

B = [0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];

% Controllability Matrix for checking whether the system is controllable or not.
Contmatrix = [B A*B A*A*B A*A*A*B A*A*A*A*B A*A*A*A*A*B];

% checking the rank of the controllability matrix, matrix should be full rank i.e 6
rank_cont = rank(Contmatrix);
disp("Rank:");
disp(rank_cont);

if (rank(Contmatrix)==6)
    disp("Rank of Controllability matrix matches the order of A, so the System is controllable")
else
    disp("Rank of Controllability matrix does not match the order of A, so the System is uncontrollable")
end

%Respnse to initial conditions when applied to linearized system.

% The initial conditions are as follows.
X_initial = [0; 0; 6; 0; 9; 0];

% Set the duration to 60 seconds
duration = 60;

% Assumptions: C matrix is represented as the output matrix, which will
% make D=0

C = eye(6); D = 0;
sys1 = ss(A, B, C, D);

% inbuilt function in MATLAB to check the initial response of the system.
figure
initial(sys1, X_initial, duration)
grid on


disp("Response of the linearized system when an LQR controller is obtained:")
% We assume the values of Q and R.
Q = [1500 0 0 0 0 0;
     0 1500 0 0 0 0;
     0 0 1500 0 0 0;
     0 0 0 1500 0 0;
     0 0 0 0 1500 0;
     0 0 0 0 0 1500];
R = 0.001;

% Q and R are the initial cost function defined in the LQR controller
% We use both the cost functions to get over the tradeoff
% inbuilt MATLAB function for designing an LQR controller
[K_Gain_matrix, Po_def_Ric, Poles] = lqr(A,B,Q,R);


sys_2 = ss(A-(B*K_Gain_matrix),B,C,D);
% Using the K matrix to define ss
figure
initial(sys_2,X_initial)
grid on


% The tweaked conditions are as follows.
X_tweaked = [0; 0; 15; 0; 25; 0];
% We assume the values of Q and R.
Q2 = [250 0 0 0 0 0;
     0 250 0 0 0 0;
     0 0 250 0 0 0;
     0 0 0 250 0 0;
     0 0 0 0 250 0;
     0 0 0 0 0 250];
R2 = 0.01;

% Q and R are the tweaked cost functions defined in the LQR controller
% We use both cost functions to get over the tradeoff
% inbuilt MATLAB function for designing an LQR controller
[K_Gain_matrix2, Po_def_Ric2, Poles2] = lqr(A,B,Q2,R2);


sys_3 = ss(A-(B*K_Gain_matrix2),B,C,D);
% Using the K matrix to define ss
figure
initial(sys_3,X_tweaked)
grid on





