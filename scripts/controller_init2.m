%#ok<*UNRCH> 

% trapezoidal, cosine, cubic
MP = 'trapezoidal';

% len = 30;
% tm = 8;
% k = 0.3552;
%
% max_vel = 8.9;
% max_acc = 10;

len = 30;
tm = 5;
k = 0.45;

max_vel = 10;
max_acc = 8;
max_dec = 5;

t1 = 0;
t3 = 0;

ta = 0;


% --------------------------- ONLY ONE SETTING SHOULD BE ON
USE_ACC_K = 0;
USE_ACC_TM = 0;
USE_ACC = 0;

USE_VEL_K = 0;
USE_VEL_TM = 0;
USE_VEL = 0;

USE_ACC_VEL = 0;

USE_K = 0;
USE_TM = 0;

USE_K_TM = 1;
% --------------------------- ONLY ONE SETTING SHOULD BE ON


% (6*D^2*R*tm - 6*D^2*R*ta - 3*A1*D*Ke^2*t1^2 - 3*A2*D*Ke^2*t3^2 - 6*A1*D*Ke^2*ta^2 + 6*A^2*A1^2*R*t1 + 6*A^2*A2^2*R*t3 - 12*A^2*A1^2*R*ta + 3*A*A1^2*Ke^2*t1^2 + 3*A*A2^2*Ke^2*t3^2 - 4*A1^2*B*Ke^2*t1^3 + 2*A2^2*B*Ke^2*t3^3 + 2*A1^2*B*Ke^2*ta^3 + 3*A1^2*B^2*L*t1^2 + 3*A2^2*B^2*L*t3^2 - 4*A1^2*B^2*R*t1^3 + 2*A2^2*B^2*R*t3^3 + 2*A1^2*B^2*R*ta^3 - 12*A1^2*B*Ke^2*t1*ta^2 + 12*A1^2*B*Ke^2*t1^2*ta + 6*A1^2*B*Ke^2*t1^2*tm + 6*A1^2*B*Ke^2*ta^2*tm - 12*A1^2*B^2*R*t1*ta^2 + 12*A1^2*B^2*R*t1^2*ta + 6*A1^2*B^2*R*t1^2*tm + 6*A1^2*B^2*R*ta^2*tm + 6*A*A1^2*B*L*t1 + 6*A*A2^2*B*L*t3 - 12*A*A1^2*B*L*ta - 6*A1*B*D*R*t1^2 - 6*A2*B*D*R*t3^2 - 12*A1*B*D*R*ta^2 + 6*A1*D*Ke^2*t1*ta + 6*A1*D*Ke^2*t1*tm - 6*A1*D*Ke^2*ta*tm + 6*A*A1^2*B*R*t1^2 + 6*A*A2^2*B*R*t3^2 - 6*A*A1^2*Ke^2*t1*ta - 6*A1^2*B^2*L*t1*ta + 6*A1*B*D*L*t1 - 6*A2*B*D*L*t3 - 12*A1*B*D*L*ta + 12*A*A1*D*R*t1 - 12*A*A2*D*R*t3 - 24*A*A1*D*R*ta - 6*A1*A2*B*Ke^2*t1*t3^2 + 6*A1*A2*B*Ke^2*t3^2*ta - 6*A1*A2*B^2*R*t1*t3^2 + 6*A1*A2*B^2*R*t3^2*ta - 12*A1^2*B*Ke^2*t1*ta*tm - 12*A1^2*B^2*R*t1*ta*tm + 12*A1*B*D*R*t1*ta + 12*A1*B*D*R*t1*tm - 12*A1*B*D*R*ta*tm - 6*A*A1*A2*Ke^2*t1*t3 + 6*A*A1*A2*Ke^2*t3*ta - 6*A1*A2*B^2*L*t1*t3 + 6*A1*A2*B^2*L*t3*ta - 12*A*A1^2*B*R*t1*ta - 12*A*A1*A2*B*R*t1*t3 + 12*A*A1*A2*B*R*t3*ta)/(6*Ke^2)
% -(2*A^2*L*Vm^2*t1^2*pi^3 - 2*A^2*R*Vm^2*t1^2*pi^3 + 8*pi*D*Ke^2*Vm*t1^3*t3 - 8*D*Ke^2*Vm*t1^3*t3*sin((pi*ta)/t1) + 16*pi*D^2*R*t1^2*t3*ta - 16*pi*D^2*R*t1^2*t3*tm + 8*D*Ke^2*Vm*t1^3*t3*sin((pi*(t1 - ta))/t1) + 8*pi*A*Ke^2*Vm^2*t1^2*t3 + 2*A^2*L*Vm^2*t1*t3*pi^3 - 4*A^2*L*Vm^2*t3*ta*pi^3 + 10*pi*B*Ke^2*Vm^2*t1^3*t3 + 8*pi*B^2*L*Vm^2*t1^2*t3 + 8*pi*D*Ke^2*Vm*t1^2*t3^2 - 2*A^2*R*Vm^2*t1*t3*pi^3 + 4*A^2*R*Vm^2*t3*ta*pi^3 + 10*pi*B^2*R*Vm^2*t1^3*t3 - 8*B*Ke^2*Vm^2*t1^3*t3*sin((pi*ta)/t1) + B*Ke^2*Vm^2*t1^3*t3*sin((2*pi*ta)/t1) - 8*B^2*R*Vm^2*t1^3*t3*sin((pi*ta)/t1) + B^2*R*Vm^2*t1^3*t3*sin((2*pi*ta)/t1) + 8*B*Ke^2*Vm^2*t1^3*t3*sin((pi*(t1 - ta))/t1) - B*Ke^2*Vm^2*t1^3*t3*sin((2*pi*(t1 - ta))/t1) + 8*B^2*R*Vm^2*t1^3*t3*sin((pi*(t1 - ta))/t1) - B^2*R*Vm^2*t1^3*t3*sin((2*pi*(t1 - ta))/t1) + 10*pi*B*Ke^2*Vm^2*t1^2*t3^2 + 10*pi*B^2*R*Vm^2*t1^2*t3^2 + 4*A^2*L*Vm^2*t1*t3*pi^2*sin((pi*ta)/t1) - A^2*L*Vm^2*t1*t3*pi^2*sin((2*pi*ta)/t1) - A^2*R*Vm^2*t1*t3*pi^2*sin((2*pi*ta)/t1) - 4*pi*B*Ke^2*Vm^2*t1^2*t3*ta - 16*pi*B*Ke^2*Vm^2*t1^2*t3*tm - 4*pi*B^2*R*Vm^2*t1^2*t3*ta - 16*pi*B^2*R*Vm^2*t1^2*t3*tm + 16*pi*B*D*L*Vm*t1^2*t3 + 32*pi*A*D*R*Vm*t1^2*t3 + 16*pi*B*D*R*Vm*t1^3*t3 + 4*A*Ke^2*Vm^2*t1^2*t3*pi*cos((pi*(t1 - ta))/t1) - A*Ke^2*Vm^2*t1^2*t3*pi*cos((2*pi*(t1 - ta))/t1) + 4*B^2*L*Vm^2*t1^2*t3*pi*cos((pi*(t1 - ta))/t1) - B^2*L*Vm^2*t1^2*t3*pi*cos((2*pi*(t1 - ta))/t1) - 4*A^2*L*Vm^2*t1*t3*pi^2*sin((pi*(t1 - ta))/t1) + A^2*L*Vm^2*t1*t3*pi^2*sin((2*pi*(t1 - ta))/t1) + A^2*R*Vm^2*t1*t3*pi^2*sin((2*pi*(t1 - ta))/t1) - 16*B*D*R*Vm*t1^3*t3*sin((pi*ta)/t1) + 16*B*D*R*Vm*t1^3*t3*sin((pi*(t1 - ta))/t1) + 16*pi*A*B*R*Vm^2*t1^2*t3 + 16*pi*B*D*R*Vm*t1^2*t3^2 - 16*pi*D*Ke^2*Vm*t1^2*t3*tm - 4*A*Ke^2*Vm^2*t1^2*t3*pi*cos((pi*ta)/t1) + A*Ke^2*Vm^2*t1^2*t3*pi*cos((2*pi*ta)/t1) - 4*B^2*L*Vm^2*t1^2*t3*pi*cos((pi*ta)/t1) + B^2*L*Vm^2*t1^2*t3*pi*cos((2*pi*ta)/t1) + 8*A*B*R*Vm^2*t1^2*t3*pi*cos((pi*(t1 - ta))/t1) - 2*A*B*R*Vm^2*t1^2*t3*pi*cos((2*pi*(t1 - ta))/t1) - 4*A*B*L*Vm^2*t1*t3*pi^2*sin((pi*(t1 - ta))/t1) + 2*A*B*L*Vm^2*t1*t3*pi^2*sin((2*pi*(t1 - ta))/t1) - 8*B*D*L*Vm*t1^2*t3*pi*cos((pi*ta)/t1) - 16*A*D*R*Vm*t1^2*t3*pi*cos((pi*ta)/t1) + 8*A*D*L*Vm*t1*t3*pi^2*sin((pi*ta)/t1) - 32*pi*B*D*R*Vm*t1^2*t3*tm + 8*B*D*L*Vm*t1^2*t3*pi*cos((pi*(t1 - ta))/t1) + 16*A*D*R*Vm*t1^2*t3*pi*cos((pi*(t1 - ta))/t1) - 8*A*D*L*Vm*t1*t3*pi^2*sin((pi*(t1 - ta))/t1) - 8*A*B*R*Vm^2*t1^2*t3*pi*cos((pi*ta)/t1) + 2*A*B*R*Vm^2*t1^2*t3*pi*cos((2*pi*ta)/t1) + 4*A*B*L*Vm^2*t1*t3*pi^2*sin((pi*ta)/t1) - 2*A*B*L*Vm^2*t1*t3*pi^2*sin((2*pi*ta)/t1))/(16*Ke^2*t1^2*t3*pi)
% -(20*A^2*L*Vm^2*t1^4 - 20*A^2*R*Vm^2*t1^4 + 15*D^2*R*t1^4*t3*ta - 15*D^2*R*t1^4*t3*tm + 20*A^2*L*Vm^2*t1^3*t3 + 20*A^2*L*Vm^2*t3*ta^3 + 7*B*Ke^2*Vm^2*t1^5*t3 + 6*B*Ke^2*Vm^2*t3*ta^5 + 5*D*Ke^2*Vm*t1^4*t3^2 - 20*A^2*R*Vm^2*t1^3*t3 + 40*A^2*R*Vm^2*t3*ta^3 + 7*B^2*R*Vm^2*t1^5*t3 + 6*B^2*R*Vm^2*t3*ta^5 + 30*A*D*L*Vm*t1^4 + 7*B*Ke^2*Vm^2*t1^4*t3^2 + 7*B^2*R*Vm^2*t1^4*t3^2 + 5*D*Ke^2*Vm*t1^5*t3 + 45*A*Ke^2*Vm^2*t1^2*t3*ta^2 + 10*B*Ke^2*Vm^2*t1^2*t3*ta^3 + 45*B^2*L*Vm^2*t1^2*t3*ta^2 + 10*B^2*R*Vm^2*t1^2*t3*ta^3 + 60*A*B*L*Vm^2*t3*ta^3 + 10*B*D*R*Vm*t1^4*t3^2 - 15*D*Ke^2*Vm*t1^4*t3*tm - 30*A*Ke^2*Vm^2*t1*t3*ta^3 - 30*A^2*L*Vm^2*t1*t3*ta^2 - 30*A^2*L*Vm^2*t1^2*t3*ta - 15*B*Ke^2*Vm^2*t1*t3*ta^4 - 15*B*Ke^2*Vm^2*t1^4*t3*tm - 30*B^2*L*Vm^2*t1*t3*ta^3 - 10*D*Ke^2*Vm*t1^2*t3*ta^3 + 15*D*Ke^2*Vm*t1^3*t3*ta^2 - 60*A^2*R*Vm^2*t1*t3*ta^2 + 60*A^2*R*Vm^2*t1^2*t3*ta - 15*B^2*R*Vm^2*t1*t3*ta^4 - 15*B^2*R*Vm^2*t1^4*t3*tm + 30*A*D*L*Vm*t1^3*t3 + 10*B*D*R*Vm*t1^5*t3 + 90*A*B*R*Vm^2*t1^2*t3*ta^2 - 60*A*D*L*Vm*t1^2*t3*ta + 30*B*D*L*Vm*t1^3*t3*ta + 60*A*D*R*Vm*t1^3*t3*ta - 30*B*D*R*Vm*t1^4*t3*tm - 90*A*B*L*Vm^2*t1*t3*ta^2 + 30*A*B*L*Vm^2*t1^2*t3*ta - 60*A*B*R*Vm^2*t1*t3*ta^3 - 20*B*D*R*Vm*t1^2*t3*ta^3 + 30*B*D*R*Vm*t1^3*t3*ta^2)/(15*Ke^2*t1^4*t3)

% ((6*D^2*R*x(1) - 6*D^2*R*ta - 3*(x(4)/x(2))*D*Ke^2*x(2)^2 - 3*(x(4)/x(3))*D*Ke^2*x(3)^2 - 6*(x(4)/x(2))*D*Ke^2*ta^2 + 6*A^2*(x(4)/x(2))^2*R*x(2) + 6*A^2*(x(4)/x(3))^2*R*x(3) - 12*A^2*(x(4)/x(2))^2*R*ta + 3*A*(x(4)/x(2))^2*Ke^2*x(2)^2 + 3*A*(x(4)/x(3))^2*Ke^2*x(3)^2 - 4*(x(4)/x(2))^2*B*Ke^2*x(2)^3 + 2*(x(4)/x(3))^2*B*Ke^2*x(3)^3 + 2*(x(4)/x(2))^2*B*Ke^2*ta^3 + 3*(x(4)/x(2))^2*B^2*L*x(2)^2 + 3*(x(4)/x(3))^2*B^2*L*x(3)^2 - 4*(x(4)/x(2))^2*B^2*R*x(2)^3 + 2*(x(4)/x(3))^2*B^2*R*x(3)^3 + 2*(x(4)/x(2))^2*B^2*R*ta^3 - 12*(x(4)/x(2))^2*B*Ke^2*x(2)*ta^2 + 12*(x(4)/x(2))^2*B*Ke^2*x(2)^2*ta + 6*(x(4)/x(2))^2*B*Ke^2*x(2)^2*x(1) + 6*(x(4)/x(2))^2*B*Ke^2*ta^2*x(1) - 12*(x(4)/x(2))^2*B^2*R*x(2)*ta^2 + 12*(x(4)/x(2))^2*B^2*R*x(2)^2*ta + 6*(x(4)/x(2))^2*B^2*R*x(2)^2*x(1) + 6*(x(4)/x(2))^2*B^2*R*ta^2*x(1) + 6*A*(x(4)/x(2))^2*B*L*x(2) + 6*A*(x(4)/x(3))^2*B*L*x(3) - 12*A*(x(4)/x(2))^2*B*L*ta - 6*(x(4)/x(2))*B*D*R*x(2)^2 - 6*(x(4)/x(3))*B*D*R*x(3)^2 - 12*(x(4)/x(2))*B*D*R*ta^2 + 6*(x(4)/x(2))*D*Ke^2*x(2)*ta + 6*(x(4)/x(2))*D*Ke^2*x(2)*x(1) - 6*(x(4)/x(2))*D*Ke^2*ta*x(1) + 6*A*(x(4)/x(2))^2*B*R*x(2)^2 + 6*A*(x(4)/x(3))^2*B*R*x(3)^2 - 6*A*(x(4)/x(2))^2*Ke^2*x(2)*ta - 6*(x(4)/x(2))^2*B^2*L*x(2)*ta + 6*(x(4)/x(2))*B*D*L*x(2) - 6*(x(4)/x(3))*B*D*L*x(3) - 12*(x(4)/x(2))*B*D*L*ta + 12*A*(x(4)/x(2))*D*R*x(2) - 12*A*(x(4)/x(3))*D*R*x(3) - 24*A*(x(4)/x(2))*D*R*ta - 6*(x(4)/x(2))*(x(4)/x(3))*B*Ke^2*x(2)*x(3)^2 + 6*(x(4)/x(2))*(x(4)/x(3))*B*Ke^2*x(3)^2*ta - 6*(x(4)/x(2))*(x(4)/x(3))*B^2*R*x(2)*x(3)^2 + 6*(x(4)/x(2))*(x(4)/x(3))*B^2*R*x(3)^2*ta - 12*(x(4)/x(2))^2*B*Ke^2*x(2)*ta*x(1) - 12*(x(4)/x(2))^2*B^2*R*x(2)*ta*x(1) + 12*(x(4)/x(2))*B*D*R*x(2)*ta + 12*(x(4)/x(2))*B*D*R*x(2)*x(1) - 12*(x(4)/x(2))*B*D*R*ta*x(1) - 6*A*(x(4)/x(2))*(x(4)/x(3))*Ke^2*x(2)*x(3) + 6*A*(x(4)/x(2))*(x(4)/x(3))*Ke^2*x(3)*ta - 6*(x(4)/x(2))*(x(4)/x(3))*B^2*L*x(2)*x(3) + 6*(x(4)/x(2))*(x(4)/x(3))*B^2*L*x(3)*ta - 12*A*(x(4)/x(2))^2*B*R*x(2)*ta - 12*A*(x(4)/x(2))*(x(4)/x(3))*B*R*x(2)*x(3) + 12*A*(x(4)/x(2))*(x(4)/x(3))*B*R*x(3)*ta)/(6*Ke^2))
% (-(2*A^2*L*x(4)^2*x(2)^2*pi^3 - 2*A^2*R*x(4)^2*x(2)^2*pi^3 + 8*pi*D*Ke^2*x(4)*x(2)^3*x(3) - 8*D*Ke^2*x(4)*x(2)^3*x(3)*sin((pi*ta)/x(2)) + 16*pi*D^2*R*x(2)^2*x(3)*ta - 16*pi*D^2*R*x(2)^2*x(3)*x(1) + 8*D*Ke^2*x(4)*x(2)^3*x(3)*sin((pi*(x(2) - ta))/x(2)) + 8*pi*A*Ke^2*x(4)^2*x(2)^2*x(3) + 2*A^2*L*x(4)^2*x(2)*x(3)*pi^3 - 4*A^2*L*x(4)^2*x(3)*ta*pi^3 + 10*pi*B*Ke^2*x(4)^2*x(2)^3*x(3) + 8*pi*B^2*L*x(4)^2*x(2)^2*x(3) + 8*pi*D*Ke^2*x(4)*x(2)^2*x(3)^2 - 2*A^2*R*x(4)^2*x(2)*x(3)*pi^3 + 4*A^2*R*x(4)^2*x(3)*ta*pi^3 + 10*pi*B^2*R*x(4)^2*x(2)^3*x(3) - 8*B*Ke^2*x(4)^2*x(2)^3*x(3)*sin((pi*ta)/x(2)) + B*Ke^2*x(4)^2*x(2)^3*x(3)*sin((2*pi*ta)/x(2)) - 8*B^2*R*x(4)^2*x(2)^3*x(3)*sin((pi*ta)/x(2)) + B^2*R*x(4)^2*x(2)^3*x(3)*sin((2*pi*ta)/x(2)) + 8*B*Ke^2*x(4)^2*x(2)^3*x(3)*sin((pi*(x(2) - ta))/x(2)) - B*Ke^2*x(4)^2*x(2)^3*x(3)*sin((2*pi*(x(2) - ta))/x(2)) + 8*B^2*R*x(4)^2*x(2)^3*x(3)*sin((pi*(x(2) - ta))/x(2)) - B^2*R*x(4)^2*x(2)^3*x(3)*sin((2*pi*(x(2) - ta))/x(2)) + 10*pi*B*Ke^2*x(4)^2*x(2)^2*x(3)^2 + 10*pi*B^2*R*x(4)^2*x(2)^2*x(3)^2 + 4*A^2*L*x(4)^2*x(2)*x(3)*pi^2*sin((pi*ta)/x(2)) - A^2*L*x(4)^2*x(2)*x(3)*pi^2*sin((2*pi*ta)/x(2)) - A^2*R*x(4)^2*x(2)*x(3)*pi^2*sin((2*pi*ta)/x(2)) - 4*pi*B*Ke^2*x(4)^2*x(2)^2*x(3)*ta - 16*pi*B*Ke^2*x(4)^2*x(2)^2*x(3)*x(1) - 4*pi*B^2*R*x(4)^2*x(2)^2*x(3)*ta - 16*pi*B^2*R*x(4)^2*x(2)^2*x(3)*x(1) + 16*pi*B*D*L*x(4)*x(2)^2*x(3) + 32*pi*A*D*R*x(4)*x(2)^2*x(3) + 16*pi*B*D*R*x(4)*x(2)^3*x(3) + 4*A*Ke^2*x(4)^2*x(2)^2*x(3)*pi*cos((pi*(x(2) - ta))/x(2)) - A*Ke^2*x(4)^2*x(2)^2*x(3)*pi*cos((2*pi*(x(2) - ta))/x(2)) + 4*B^2*L*x(4)^2*x(2)^2*x(3)*pi*cos((pi*(x(2) - ta))/x(2)) - B^2*L*x(4)^2*x(2)^2*x(3)*pi*cos((2*pi*(x(2) - ta))/x(2)) - 4*A^2*L*x(4)^2*x(2)*x(3)*pi^2*sin((pi*(x(2) - ta))/x(2)) + A^2*L*x(4)^2*x(2)*x(3)*pi^2*sin((2*pi*(x(2) - ta))/x(2)) + A^2*R*x(4)^2*x(2)*x(3)*pi^2*sin((2*pi*(x(2) - ta))/x(2)) - 16*B*D*R*x(4)*x(2)^3*x(3)*sin((pi*ta)/x(2)) + 16*B*D*R*x(4)*x(2)^3*x(3)*sin((pi*(x(2) - ta))/x(2)) + 16*pi*A*B*R*x(4)^2*x(2)^2*x(3) + 16*pi*B*D*R*x(4)*x(2)^2*x(3)^2 - 16*pi*D*Ke^2*x(4)*x(2)^2*x(3)*x(1) - 4*A*Ke^2*x(4)^2*x(2)^2*x(3)*pi*cos((pi*ta)/x(2)) + A*Ke^2*x(4)^2*x(2)^2*x(3)*pi*cos((2*pi*ta)/x(2)) - 4*B^2*L*x(4)^2*x(2)^2*x(3)*pi*cos((pi*ta)/x(2)) + B^2*L*x(4)^2*x(2)^2*x(3)*pi*cos((2*pi*ta)/x(2)) + 8*A*B*R*x(4)^2*x(2)^2*x(3)*pi*cos((pi*(x(2) - ta))/x(2)) - 2*A*B*R*x(4)^2*x(2)^2*x(3)*pi*cos((2*pi*(x(2) - ta))/x(2)) - 4*A*B*L*x(4)^2*x(2)*x(3)*pi^2*sin((pi*(x(2) - ta))/x(2)) + 2*A*B*L*x(4)^2*x(2)*x(3)*pi^2*sin((2*pi*(x(2) - ta))/x(2)) - 8*B*D*L*x(4)*x(2)^2*x(3)*pi*cos((pi*ta)/x(2)) - 16*A*D*R*x(4)*x(2)^2*x(3)*pi*cos((pi*ta)/x(2)) + 8*A*D*L*x(4)*x(2)*x(3)*pi^2*sin((pi*ta)/x(2)) - 32*pi*B*D*R*x(4)*x(2)^2*x(3)*x(1) + 8*B*D*L*x(4)*x(2)^2*x(3)*pi*cos((pi*(x(2) - ta))/x(2)) + 16*A*D*R*x(4)*x(2)^2*x(3)*pi*cos((pi*(x(2) - ta))/x(2)) - 8*A*D*L*x(4)*x(2)*x(3)*pi^2*sin((pi*(x(2) - ta))/x(2)) - 8*A*B*R*x(4)^2*x(2)^2*x(3)*pi*cos((pi*ta)/x(2)) + 2*A*B*R*x(4)^2*x(2)^2*x(3)*pi*cos((2*pi*ta)/x(2)) + 4*A*B*L*x(4)^2*x(2)*x(3)*pi^2*sin((pi*ta)/x(2)) - 2*A*B*L*x(4)^2*x(2)*x(3)*pi^2*sin((2*pi*ta)/x(2)))/(16*Ke^2*x(2)^2*x(3)*pi))
% (-(20*A^2*L*x(4)^2*x(2)^4 - 20*A^2*R*x(4)^2*x(2)^4 + 15*D^2*R*x(2)^4*x(3)*ta - 15*D^2*R*x(2)^4*x(3)*x(1) + 20*A^2*L*x(4)^2*x(2)^3*x(3) + 20*A^2*L*x(4)^2*x(3)*ta^3 + 7*B*Ke^2*x(4)^2*x(2)^5*x(3) + 6*B*Ke^2*x(4)^2*x(3)*ta^5 + 5*D*Ke^2*x(4)*x(2)^4*x(3)^2 - 20*A^2*R*x(4)^2*x(2)^3*x(3) + 40*A^2*R*x(4)^2*x(3)*ta^3 + 7*B^2*R*x(4)^2*x(2)^5*x(3) + 6*B^2*R*x(4)^2*x(3)*ta^5 + 30*A*D*L*x(4)*x(2)^4 + 7*B*Ke^2*x(4)^2*x(2)^4*x(3)^2 + 7*B^2*R*x(4)^2*x(2)^4*x(3)^2 + 5*D*Ke^2*x(4)*x(2)^5*x(3) + 45*A*Ke^2*x(4)^2*x(2)^2*x(3)*ta^2 + 10*B*Ke^2*x(4)^2*x(2)^2*x(3)*ta^3 + 45*B^2*L*x(4)^2*x(2)^2*x(3)*ta^2 + 10*B^2*R*x(4)^2*x(2)^2*x(3)*ta^3 + 60*A*B*L*x(4)^2*x(3)*ta^3 + 10*B*D*R*x(4)*x(2)^4*x(3)^2 - 15*D*Ke^2*x(4)*x(2)^4*x(3)*x(1) - 30*A*Ke^2*x(4)^2*x(2)*x(3)*ta^3 - 30*A^2*L*x(4)^2*x(2)*x(3)*ta^2 - 30*A^2*L*x(4)^2*x(2)^2*x(3)*ta - 15*B*Ke^2*x(4)^2*x(2)*x(3)*ta^4 - 15*B*Ke^2*x(4)^2*x(2)^4*x(3)*x(1) - 30*B^2*L*x(4)^2*x(2)*x(3)*ta^3 - 10*D*Ke^2*x(4)*x(2)^2*x(3)*ta^3 + 15*D*Ke^2*x(4)*x(2)^3*x(3)*ta^2 - 60*A^2*R*x(4)^2*x(2)*x(3)*ta^2 + 60*A^2*R*x(4)^2*x(2)^2*x(3)*ta - 15*B^2*R*x(4)^2*x(2)*x(3)*ta^4 - 15*B^2*R*x(4)^2*x(2)^4*x(3)*x(1) + 30*A*D*L*x(4)*x(2)^3*x(3) + 10*B*D*R*x(4)*x(2)^5*x(3) + 90*A*B*R*x(4)^2*x(2)^2*x(3)*ta^2 - 60*A*D*L*x(4)*x(2)^2*x(3)*ta + 30*B*D*L*x(4)*x(2)^3*x(3)*ta + 60*A*D*R*x(4)*x(2)^3*x(3)*ta - 30*B*D*R*x(4)*x(2)^4*x(3)*x(1) - 90*A*B*L*x(4)^2*x(2)*x(3)*ta^2 + 30*A*B*L*x(4)^2*x(2)^2*x(3)*ta - 60*A*B*R*x(4)^2*x(2)*x(3)*ta^3 - 20*B*D*R*x(4)*x(2)^2*x(3)*ta^3 + 30*B*D*R*x(4)*x(2)^3*x(3)*ta^2)/(15*Ke^2*x(2)^4*x(3)))


if(USE_ACC_K+USE_ACC_TM+USE_ACC+USE_VEL_K+USE_VEL_TM+USE_VEL+USE_ACC_VEL+USE_K+USE_TM+USE_K_TM ~= 1)
    error('INVALID SETTING');
end

mp = 0;
if(strcmp(MP,'trapezoidal') == 1)
    mp = 1;
    funE = @(x)((6*D^2*R*x(1) - 6*D^2*R*ta - 3*(x(4)/x(2))*D*Ke^2*x(2)^2 - 3*(x(4)/x(3))*D*Ke^2*x(3)^2 - 6*(x(4)/x(2))*D*Ke^2*ta^2 + 6*A^2*(x(4)/x(2))^2*R*x(2) + 6*A^2*(x(4)/x(3))^2*R*x(3) - 12*A^2*(x(4)/x(2))^2*R*ta + 3*A*(x(4)/x(2))^2*Ke^2*x(2)^2 + 3*A*(x(4)/x(3))^2*Ke^2*x(3)^2 - 4*(x(4)/x(2))^2*B*Ke^2*x(2)^3 + 2*(x(4)/x(3))^2*B*Ke^2*x(3)^3 + 2*(x(4)/x(2))^2*B*Ke^2*ta^3 + 3*(x(4)/x(2))^2*B^2*L*x(2)^2 + 3*(x(4)/x(3))^2*B^2*L*x(3)^2 - 4*(x(4)/x(2))^2*B^2*R*x(2)^3 + 2*(x(4)/x(3))^2*B^2*R*x(3)^3 + 2*(x(4)/x(2))^2*B^2*R*ta^3 - 12*(x(4)/x(2))^2*B*Ke^2*x(2)*ta^2 + 12*(x(4)/x(2))^2*B*Ke^2*x(2)^2*ta + 6*(x(4)/x(2))^2*B*Ke^2*x(2)^2*x(1) + 6*(x(4)/x(2))^2*B*Ke^2*ta^2*x(1) - 12*(x(4)/x(2))^2*B^2*R*x(2)*ta^2 + 12*(x(4)/x(2))^2*B^2*R*x(2)^2*ta + 6*(x(4)/x(2))^2*B^2*R*x(2)^2*x(1) + 6*(x(4)/x(2))^2*B^2*R*ta^2*x(1) + 6*A*(x(4)/x(2))^2*B*L*x(2) + 6*A*(x(4)/x(3))^2*B*L*x(3) - 12*A*(x(4)/x(2))^2*B*L*ta - 6*(x(4)/x(2))*B*D*R*x(2)^2 - 6*(x(4)/x(3))*B*D*R*x(3)^2 - 12*(x(4)/x(2))*B*D*R*ta^2 + 6*(x(4)/x(2))*D*Ke^2*x(2)*ta + 6*(x(4)/x(2))*D*Ke^2*x(2)*x(1) - 6*(x(4)/x(2))*D*Ke^2*ta*x(1) + 6*A*(x(4)/x(2))^2*B*R*x(2)^2 + 6*A*(x(4)/x(3))^2*B*R*x(3)^2 - 6*A*(x(4)/x(2))^2*Ke^2*x(2)*ta - 6*(x(4)/x(2))^2*B^2*L*x(2)*ta + 6*(x(4)/x(2))*B*D*L*x(2) - 6*(x(4)/x(3))*B*D*L*x(3) - 12*(x(4)/x(2))*B*D*L*ta + 12*A*(x(4)/x(2))*D*R*x(2) - 12*A*(x(4)/x(3))*D*R*x(3) - 24*A*(x(4)/x(2))*D*R*ta - 6*(x(4)/x(2))*(x(4)/x(3))*B*Ke^2*x(2)*x(3)^2 + 6*(x(4)/x(2))*(x(4)/x(3))*B*Ke^2*x(3)^2*ta - 6*(x(4)/x(2))*(x(4)/x(3))*B^2*R*x(2)*x(3)^2 + 6*(x(4)/x(2))*(x(4)/x(3))*B^2*R*x(3)^2*ta - 12*(x(4)/x(2))^2*B*Ke^2*x(2)*ta*x(1) - 12*(x(4)/x(2))^2*B^2*R*x(2)*ta*x(1) + 12*(x(4)/x(2))*B*D*R*x(2)*ta + 12*(x(4)/x(2))*B*D*R*x(2)*x(1) - 12*(x(4)/x(2))*B*D*R*ta*x(1) - 6*A*(x(4)/x(2))*(x(4)/x(3))*Ke^2*x(2)*x(3) + 6*A*(x(4)/x(2))*(x(4)/x(3))*Ke^2*x(3)*ta - 6*(x(4)/x(2))*(x(4)/x(3))*B^2*L*x(2)*x(3) + 6*(x(4)/x(2))*(x(4)/x(3))*B^2*L*x(3)*ta - 12*A*(x(4)/x(2))^2*B*R*x(2)*ta - 12*A*(x(4)/x(2))*(x(4)/x(3))*B*R*x(2)*x(3) + 12*A*(x(4)/x(2))*(x(4)/x(3))*B*R*x(3)*ta)/(6*Ke^2));
    if(USE_ACC_K)
        c = @(x) [];
        ceq = @(x)[2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; (x(1)-x(2)-x(3))/x(1) - k; x(2)+x(3)+k*x(1)-x(1); x(4)/max_acc-x(2); x(4)/max_dec-x(3);];
        nonlinfcn = @(x)deal(c(x),ceq(x));
        lb = [0.1;0.1;0.1;0.1];
        ub = [50;50;50;50];
    elseif(USE_ACC_TM) 
        c = @(x) x(2)+x(3)-x(1);
%          max_acc*x(2)*(2*x(1)-x(2)-x(3))/2-len; max_dec*x(3)*(2*x(1)-x(2)-x(3))/2-len;
        ceq = @(x) [2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1;x(4)/max_acc-x(2); x(4)/max_dec-x(3);];
        nonlinfcn = @(x)deal(c(x),ceq(x));
        lb = [tm;0.1;0.1;0.1];
        ub = [tm;50;50;50];
    elseif(USE_ACC)
        c = @(x) x(2)+x(3)-x(1);
        ceq = @(x)[2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1;x(4)/max_acc-x(2); x(4)/max_dec-x(3);];
        nonlinfcn = @(x)deal(c(x),ceq(x));
        lb = [0.1;0.1;0.1;0.1];
        ub = [50;50;50;50];
    elseif(USE_VEL_K)
        c = @(x) x(2)+x(3)-x(1);
        ceq = @(x)[2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; (x(1)-x(2)-x(3))/x(1) - k; (2/(1+k))*(len/x(1))-x(4)];
        nonlinfcn = @(x)deal(c(x),ceq(x));
        lb = [0.1;0.1;0.1;max_vel];
        ub = [50;50;50;max_vel];
    elseif(USE_VEL_TM)
        % alright
%         k = (2*len/(max_vel*tm))-1  % (FIXED)Wrong should not take abs, the value should be positive below some absolute max vel....trapz eqns allow for only same accel and decel at opt
%         if(k<0)
%             error('MAX_VEL too high to meet the motion profile requirments');
%         elseif(k>=1)
%             error('MAX_VEL too low to meet the motion profile requirements');
%         end
        if(max_vel < len/tm)
            error('MAX_VEL too low');
        end
        c = @(x) x(2)+x(3)-x(1);
        ceq = @(x) [2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1;len-x(4)*((2*x(1)-x(2)-x(3))/2)];
        nonlinfcn = @(x)deal(c(x),ceq(x));
        lb = [tm;0.1;0.1;max_vel];
        ub = [tm;tm;tm;max_vel];
    elseif(USE_VEL)
        c = @(x) x(2)+x(3)-x(1);
        ceq = @(x) 2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1;
        nonlinfcn = @(x)deal(c(x),ceq(x));
        lb = [0.1;0.1;0.1;max_vel];
        ub = [50;50;50;max_vel];
    elseif(USE_ACC_VEL)
        c = @(x) x(2)+x(3)-x(1);
        ceq = @(x) [2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; x(4)/max_acc-x(2); x(4)/max_dec-x(3)];
        nonlinfcn = @(x)deal(c(x),ceq(x));
        lb = [0.1;0.1;0.1;max_vel];
        ub = [50;50;50;max_vel];
    elseif(USE_K)
        c = @(x) x(2)+x(3)-x(1);
        ceq = @(x)[2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; (x(1)-x(2)-x(3))/x(1) - k; x(2)+x(3)+k*x(1)-x(1); (2/(1+k))*(len/x(1))-x(4)];
        nonlinfcn = @(x)deal(c(x),ceq(x));
        lb = [0.1;0.1;0.1;0.1];
        ub = [50;50;50;50];
    elseif(USE_TM)
        c = @(x) x(2)+x(3)-x(1);
        ceq = @(x)2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1;
        nonlinfcn = @(x)deal(c(x),ceq(x));
        lb = [tm;0.1;0.1;0.1];
        ub = [tm;50;50;50];
    elseif(USE_K_TM)
        c = @(x) x(2)+x(3)-x(1);
        ceq = @(x) [2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; (x(1)-x(2)-x(3))/x(1) - k; x(2)+x(3)+k*x(1)-x(1);];
        nonlinfcn = @(x)deal(c(x),ceq(x));
        lb = [tm;0.1;0.1;0.1];
        ub = [tm;50;50;50];
    end

    hybridopts = optimoptions('fmincon','Display','iter','Algorithm','interior-point','PlotFcn','optimplotfval',EnableFeasibilityMode=true,ConstraintTolerance=1e-7, OptimalityTolerance=1e-8);
    options = optimoptions('ga','PlotFcn',"gaplotbestf",'HybridFcn',{@fmincon,hybridopts},FitnessScalingFcn="fitscalingtop",SelectionFcn="selectionstochunif",MutationFcn="mutationadaptfeasible",CrossoverFcn="crossoverlaplace",MaxGenerations=800, ConstraintTolerance=1e-11, FunctionTolerance=1e-15,Display="iter", PopulationSize=100);
    [optxs,fval] = ga(funE, 4,[],[],[],[],lb,ub, nonlinfcn,options)
    tm = optxs(1); t1 = optxs(2); t3 = optxs(3); max_vel = optxs(4);
    k = (tm-t1-t3)/tm;

elseif(strcmp(MP,'cosine') == 1)
    mp = 2;
    if(USE_ACC_VEL)
        cof = (2*pi)/(max_acc*4*len);
        k = roots([(cof*max_vel^2+1) (2*cof*max_vel^2) (cof*max_vel^2-1)]);
        k = k(k<1&k>=0);
        if(isempty(k))
            error('Unsuitable MAX_VEL and MAX_ACC values');
        end
        tm = 2*len/(max_vel*(1+k));
    elseif(USE_K_TM)
        % normal case
        max_acc = (2*pi/(1-k^2))*(len/tm^2);
        max_vel = (2/(1+k)) * (len/tm);
    end
elseif(strcmp(MP,'cubic') == 1)
    mp = 3;
    if(USE_ACC_VEL)
        cof = (-12*max_vel^2)/(max_acc*9*len);
        k = roots([(cof-1) (4*cof-1) (4*cof+2)]);
        k = k(k<1&k>=0);
        if(isempty(k))
            error('Unsuitable MAX_VEL and MAX_ACC values');
        end
        tm = 3*len/(max_vel*(2+k));
    elseif(USE_K_TM)
        % normal case
        max_acc = (-12/(k^2+k-2))*(len/tm^2);
        max_vel = (3/(2+k)) * (len/tm);
    end
end


ctrl_max_acc = max_acc;
ctrl_max_dec = max_dec;
ctrl_max_vel = max_vel;

ctrl_len = len;
ctrl_tm = tm;
ctrl_k = k;

ctrl_t1 = t1;
ctrl_t3 = t3;

ctrl_mp = mp;



