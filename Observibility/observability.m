clear variables;
clc;

syms M m1 m2 l1 l2 g;

A = [0 1 0 0 0 0;
    0 0 -(g*m1)/M 0 -(g*m2)/M 0;
    0 0 0 1 0 0;
    0 0 (-g*M -m1*g)/(M*l1) 0 -(g*m2)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -(g*m1)/(M*l2) 0 (-g*M -g*m2)/(M*l2) 0];

B = [0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];

C1 = [1 0 0 0 0 0]; %Formulating with respect to x component
C2 = [0 0 1 0 0 0; 0 0 0 0 1 0]; %Formulating with respect to theta1 and theta2
C3 = [1 0 0 0 0 0; 0 0 0 0 1 0]; %Formulating with respect to x and theta2
C4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0]; %Formulating with respect to x,theta1 and theta2

Obs_check1 = rank([C1' A'*C1' ((A')^2)*C1' ((A')^3)*C1' ((A')^4)*C1' ((A')^5)*C1']);
Obs_check2 = rank([C1' A'*C2' ((A')^2)*C2' ((A')^3)*C2' ((A')^4)*C2' ((A')^5)*C2']);
Obs_check3 = rank([C1' A'*C3' ((A')^2)*C3' ((A')^3)*C3' ((A')^4)*C3' ((A')^5)*C3']);
Obs_check4 = rank([C1' A'*C4' ((A')^2)*C4' ((A')^3)*C4' ((A')^4)*C4' ((A')^5)*C4']);

disp(Obs_check1)
disp(Obs_check2)
disp(Obs_check3)
disp(Obs_check4)




