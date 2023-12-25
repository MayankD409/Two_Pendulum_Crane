clear variables;
clc;

% Defining variables

m1 = 100;
m2 = 100;
M = 1000;
l1 = 20;
l2 = 10;
g = 9.81;
t_span = 0:0.1:100;
% q = [x dx t1 dt1 t2 dt2];
%Enter initial conditions
q0 = [2 0 deg2rad(0) 0 deg2rad(0) 0];

% Linearized Model
A = [0 1 0 0 0 0;
    0 0 -(g*m1)/M 0 -(g*m2)/M 0;
    0 0 0 1 0 0;
    0 0 (-g*M -m1*g)/(M*l1) 0 -(g*m2)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -(g*m1)/(M*l2) 0 (-g*M -g*m2)/(M*l2) 0];

B = [0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];

c1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
d = [1;0;0];
sys1 = ss(A,B,c1,d);

% LQR Controller
Q = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
R = 0.1;
[K,S,P] = lqr(A,B,Q,R);
sys = ss(A-B*K,B,c1,d);
% step(sys,200);

% Kalman Estimator Design
Bd = 0.01*eye(6);                %Process Noise
Vn = 0.001;                      %Measurement Noise
[L1,P,E] = lqe(A,Bd,c1,Bd,Vn*eye(3)); %Considering vector output: x(t)
Ac1 = A-(L1*c1);
e_sys1 = ss(Ac1,[B L1],c1,0);

% Non-linear Model LQG Response
[t,q1] = ode45(@(t,q)nonLinearObs1(t,q,-K*q,L1),t_span,q0);
figure();
hold on
plot(t,q1(:,1))
ylabel('state variable')
xlabel('t')
title('Original System LQG fo x(t)')
legend('x')
hold off