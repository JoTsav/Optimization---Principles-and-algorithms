
% @author Joseph

% @note Calls \ref runCG conjugateGradient.m

Q = [4, 1, 0, 0; 1, 3, 1, 0; 0, 1, 2, 1; 0, 0, 1, 3];
b = [1; 2; 3; 4]; x0 = [0; 0; 0; 0];
[D,solution] = conjugateGradient(Q,b,x0, 1); % solution
