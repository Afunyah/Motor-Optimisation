X1 = R*D^2/Ke^2;
X2 = R*2*B*D/Ke^2 + D;
X3 = R*B^2/Ke^2 + B;
X4 = (R*A^2 + A*B*L)/Ke^2;

len = 30;
% tm = 5;
% k = 0.3334;


funE = @(x) (X1*x(1)  +  X2*len  +  (2*len/(x(1)+x(1)*x(2)))^2  * ( X3*(x(1) + 2*x(1)*x(2))/3 + X4*4/(x(1)-x(1)*x(2))));
% rng default
lb = [0 0.01];
ub = [30 0.99];
optxy = ga(funE, 2, [],[],[],[],lb, ub)