
X1 = R*D^2/Ke^2;
X2 = R*2*B*D/Ke^2 + D;
X3 = R*B^2/Ke^2 + B;
X4 = (R*A^2 + A*B*L)/Ke^2;

len = 30;
% tm = 5;
% k = 0.3334;

tm = optimvar('tm');
k = optimvar('k');

E = X1*tm  +  X2*len  +  (2*len/(tm+tm*k))^2  * ( X3*(tm + 2*tm*k)/3 + X4*4/(tm-tm*k));

prob = optimproblem("Objective",E);
x0.tm = 3;
x0.k = 0.5;
[sol,fval] = solve(prob,x0) 