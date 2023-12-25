clear variables;
clc;
% Defining the variables 
syms M m1 m2 l1 l2 g;

% state-space representation  of the linearized system
A = [0 1 0 0 0 0;
    0 0 -(g*m1)/M 0 -(g*m2)/M 0;
    0 0 0 1 0 0;
    0 0 (-g*M -m1*g)/(M*l1) 0 -(g*m2)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -(g*m1)/(M*l2) 0 (-g*M -g*m2)/(M*l2) 0];
disp("A")
disp(A)

B = [0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];

disp("B")
disp(B)

% Controllability Matrix for checking whether the system is controllable or
% not.

Contmatrix = [B A*B A*A*B A*A*A*B A*A*A*A*B A*A*A*A*A*B];
disp(Contmatrix)

% checking the rank of the controllability of the matrix, matrix should be full rank i.e 6

rank_cont = rank(Contmatrix);
disp("Rank:");
disp(rank_cont);

% Condition of controllability (to check the determinant of the matrix

% Check if the system is controllable
isControllable = rank(Contmatrix) == 6; 

% Display the result
if isControllable
    disp('The system is controllable.');
else
    disp('The system is not controllable.');
end




