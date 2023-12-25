clear variables;
clc;
% Defining variables

m1 = 100;
m2 = 100;
M = 1000;
l1 = 20;
l2 = 10;
g = 9.81;
q0 = [2 0 deg2rad(17) 0 deg2rad(30) 0];
t_span = 0:0.1:100;

% Observability Check
A = [0 1 0 0 0 0;
     0 0 -(g*m1)/M 0 -(g*m2)/M 0;
     0 0 0 1 0 0;
     0 0 (-g*M -m1*g)/(M*l1) 0 -(g*m2)/(M*l1) 0;
     0 0 0 0 0 1;
     0 0 -(g*m1)/(M*l2) 0 (-g*M -g*m2)/(M*l2) 0];

B = [0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];
c1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
c2 = [0 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
c3 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0];
c4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
d = [0; 0; 0];
Obs1 = rank([c1' A'*c1' ((A')^2)*c1' ((A')^3)*c1' ((A')^4)*c1' ((A')^5)*c1']);
Obs2 = rank([c2' A'*c2' ((A')^2)*c2' ((A')^3)*c2' ((A')^4)*c2' ((A')^5)*c2']);
Obs3 = rank([c3' A'*c3' ((A')^2)*c3' ((A')^3)*c3' ((A')^4)*c3' ((A')^5)*c3']);
Obs4 = rank([c4' A'*c4' ((A')^2)*c4' ((A')^3)*c4' ((A')^4)*c4' ((A')^5)*c4']);

sys1 = ss(A,B,c1,d);
sys3 = ss(A,B,c3,d);
sys4 = ss(A,B,c4,d);

% Kalman Estimator Design
Bd = 0.01*eye(6);                %Process Noise
Vn = 0.001;                      %Measurement Noise
[L1,P,E] = lqe(A,Bd,c1,Bd,Vn*eye(3));
[L3,P,E] = lqe(A,Bd,c3,Bd,Vn*eye(3));
[L4,P,E] = lqe(A,Bd,c4,Bd,Vn*eye(3));

Ac1 = A-(L1*c1);
Ac3 = A-(L3*c3);
Ac4 = A-(L4*c4);
e_sys1 = ss(Ac1,[B L1],c1,0);
e_sys3 = ss(Ac3,[B L3],c3,0);
e_sys4 = ss(Ac4,[B L4],c4,0);


% Generating plot for step input
u_Step = 0*t_span;
u_Step(200:length(t_span)) = 1;

[y1,t] = lsim(sys1,u_Step,t_span);
[x1,t] = lsim(e_sys1,[u_Step;y1'],t_span);

[y3,t] = lsim(sys3,u_Step,t_span);
[x3,t] = lsim(e_sys3,[u_Step;y3'],t_span);

[y4,t] = lsim(sys4,u_Step,t_span);
[x4,t] = lsim(e_sys4,[u_Step;y4'],t_span);

figure();
hold on
plot(t,y1(:,1),'r','Linewidth',2)
plot(t,x1(:,1),'k--','Linewidth',1)
ylabel('State Variables')
xlabel('t')
legend('x(t)','Estimated x(t)')
title('(x(t)')
hold off

figure();
hold on
plot(t,y3(:,1),'r','Linewidth',2)
plot(t,y3(:,3),'b','Linewidth',2)
plot(t,x3(:,1),'k--','Linewidth',1)
plot(t,x3(:,3),'m--','Linewidth',1)
ylabel('State Variables')
xlabel('t')
legend('x(t)','theta_2(t)','Estimated x(t)','Estimated theta_2(t)')
title('(x(t),theta_2(t))')
hold off

figure();
hold on
plot(t,y4(:,1),'r','Linewidth',2)
plot(t,y4(:,2),'g','Linewidth',2)
plot(t,y4(:,3),'b','Linewidth',2)
plot(t,x4(:,1),'k--','Linewidth',1)
plot(t,x4(:,2),'r--','Linewidth',1)
plot(t,x4(:,3),'m--','Linewidth',1)
ylabel('State Variables')
xlabel('t')
legend('x(t)','theta_1(t)','theta_2(t)','Estimated x(t)','Estimated theta_1(t)','Estimated theta_2(t)')
title('(x(t),theta_1(t),theta_2(t))')
hold off

% Leuenberger Observer Response for linearized system
[t,q1] = ode45(@(t,q)linearObs1(t,q,L1),t_span,q0);
figure();
hold on
plot(t,q1(:,1))
ylabel('state variables')
xlabel('t')
title('For linearized system: x(t)')
legend('x')
hold off

[t,q3] = ode45(@(t,q)linearObs3(t,q,L3),t_span,q0);
figure();
hold on
plot(t,q3(:,1))
plot(t,q3(:,5))
ylabel('state variables')
xlabel('t')
title('For linearized system: (x(t),theta_2(t))')
legend('x','theta_2')
hold off

[t,q4] = ode45(@(t,q)linearObs4(t,q,L4),t_span,q0);
figure();
hold on
plot(t,q4(:,1))
plot(t,q4(:,3))
plot(t,q4(:,5))
ylabel('state variables')
xlabel('t')
title('For linearized system: (x(t),theta_1(t),theta_2(t))')
legend('x','theta_1','theta_2')
hold off

% Leuenberger Observer Response for the original non-linear system
[t,q1] = ode45(@(t,q)nonLinearObs1(t,q,1,L1),t_span,q0);
figure();
hold on
plot(t,q1(:,1))
ylabel('state variables')
xlabel('t')
title('For Non-linear system: x(t)')
legend('x')
hold off

[t,q3] = ode45(@(t,q)nonLinearObs3(t,q,1,L3),t_span,q0);
figure();
hold on
plot(t,q3(:,1))
plot(t,q3(:,5))
ylabel('state variables')
xlabel('t')
title('For Non-linear system: (x(t),theta_2(t))')
legend('x','theta_2')
hold off

[t,q4] = ode45(@(t,q)nonLinearObs4(t,q,1,L4),t_span,q0);
figure();
hold on
plot(t,q4(:,1))
plot(t,q4(:,3))
plot(t,q4(:,5))
ylabel('state variables')
xlabel('t')
title('For Non-linear system: (x(t),theta_1(t),theta_2(t))')
legend('x','theta_1','theta_2')
hold off