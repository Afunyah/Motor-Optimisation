
% Cubic

% 1:tm 2:t1 3:t3 4:Vm

len = 30;
tm = 9;
k = 0.34;

ta = 0;

% -(x(1)-x(2)-x(3))/x(1);
%         (x(1)-x(2)-x(3))/x(1)-0.99;
%         -(3*len/(x(4)*x(1)))+2;
%         (3*len/(x(4)*x(1)))-2-0.99;

c = @(x)[
        x(2)+x(3)-x(1)
        ];
%         x(2)+x(3)-x(1)
ceq = @(x)[2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1;];
nonlinfcn = @(x)deal(c(x),ceq(x));



% funE = @(x) (-(20*A^2*L*x(4)^2*x(2)^4 - 20*A^2*R*x(4)^2*x(2)^4 + 15*D^2*R*x(2)^4*x(3)*ta - 15*D^2*R*x(2)^4*x(3)*x(1) + 20*A^2*L*x(4)^2*x(2)^3*x(3) + 20*A^2*L*x(4)^2*x(3)*ta^3 + 7*B*Ke^2*x(4)^2*x(2)^5*x(3) + 6*B*Ke^2*x(4)^2*x(3)*ta^5 + 5*D*Ke^2*x(4)*x(2)^4*x(3)^2 - 20*A^2*R*x(4)^2*x(2)^3*x(3) + 40*A^2*R*x(4)^2*x(3)*ta^3 + 7*B^2*R*x(4)^2*x(2)^5*x(3) + 6*B^2*R*x(4)^2*x(3)*ta^5 + 30*A*D*L*x(4)*x(2)^4 + 7*B*Ke^2*x(4)^2*x(2)^4*x(3)^2 + 7*B^2*R*x(4)^2*x(2)^4*x(3)^2 + 5*D*Ke^2*x(4)*x(2)^5*x(3) + 45*A*Ke^2*x(4)^2*x(2)^2*x(3)*ta^2 + 10*B*Ke^2*x(4)^2*x(2)^2*x(3)*ta^3 + 45*B^2*L*x(4)^2*x(2)^2*x(3)*ta^2 + 10*B^2*R*x(4)^2*x(2)^2*x(3)*ta^3 + 60*A*B*L*x(4)^2*x(3)*ta^3 + 10*B*D*R*x(4)*x(2)^4*x(3)^2 - 15*D*Ke^2*x(4)*x(2)^4*x(3)*x(1) - 30*A*Ke^2*x(4)^2*x(2)*x(3)*ta^3 - 30*A^2*L*x(4)^2*x(2)*x(3)*ta^2 - 30*A^2*L*x(4)^2*x(2)^2*x(3)*ta - 15*B*Ke^2*x(4)^2*x(2)*x(3)*ta^4 - 15*B*Ke^2*x(4)^2*x(2)^4*x(3)*x(1) - 30*B^2*L*x(4)^2*x(2)*x(3)*ta^3 - 10*D*Ke^2*x(4)*x(2)^2*x(3)*ta^3 + 15*D*Ke^2*x(4)*x(2)^3*x(3)*ta^2 - 60*A^2*R*x(4)^2*x(2)*x(3)*ta^2 + 60*A^2*R*x(4)^2*x(2)^2*x(3)*ta - 15*B^2*R*x(4)^2*x(2)*x(3)*ta^4 - 15*B^2*R*x(4)^2*x(2)^4*x(3)*x(1) + 30*A*D*L*x(4)*x(2)^3*x(3) + 10*B*D*R*x(4)*x(2)^5*x(3) + 90*A*B*R*x(4)^2*x(2)^2*x(3)*ta^2 - 60*A*D*L*x(4)*x(2)^2*x(3)*ta + 30*B*D*L*x(4)*x(2)^3*x(3)*ta + 60*A*D*R*x(4)*x(2)^3*x(3)*ta - 30*B*D*R*x(4)*x(2)^4*x(3)*x(1) - 90*A*B*L*x(4)^2*x(2)*x(3)*ta^2 + 30*A*B*L*x(4)^2*x(2)^2*x(3)*ta - 60*A*B*R*x(4)^2*x(2)*x(3)*ta^3 - 20*B*D*R*x(4)*x(2)^2*x(3)*ta^3 + 30*B*D*R*x(4)*x(2)^3*x(3)*ta^2)/(15*Ke^2*x(2)^4*x(3)));

funE = @(x)((6*D^2*R*x(1) - 6*D^2*R*ta - 3*(x(4)/x(2))*D*Ke^2*x(2)^2 - 3*(x(4)/x(3))*D*Ke^2*x(3)^2 - 6*(x(4)/x(2))*D*Ke^2*ta^2 + 6*A^2*(x(4)/x(2))^2*R*x(2) + 6*A^2*(x(4)/x(3))^2*R*x(3) - 12*A^2*(x(4)/x(2))^2*R*ta + 3*A*(x(4)/x(2))^2*Ke^2*x(2)^2 + 3*A*(x(4)/x(3))^2*Ke^2*x(3)^2 - 4*(x(4)/x(2))^2*B*Ke^2*x(2)^3 + 2*(x(4)/x(3))^2*B*Ke^2*x(3)^3 + 2*(x(4)/x(2))^2*B*Ke^2*ta^3 + 3*(x(4)/x(2))^2*B^2*L*x(2)^2 + 3*(x(4)/x(3))^2*B^2*L*x(3)^2 - 4*(x(4)/x(2))^2*B^2*R*x(2)^3 + 2*(x(4)/x(3))^2*B^2*R*x(3)^3 + 2*(x(4)/x(2))^2*B^2*R*ta^3 - 12*(x(4)/x(2))^2*B*Ke^2*x(2)*ta^2 + 12*(x(4)/x(2))^2*B*Ke^2*x(2)^2*ta + 6*(x(4)/x(2))^2*B*Ke^2*x(2)^2*x(1) + 6*(x(4)/x(2))^2*B*Ke^2*ta^2*x(1) - 12*(x(4)/x(2))^2*B^2*R*x(2)*ta^2 + 12*(x(4)/x(2))^2*B^2*R*x(2)^2*ta + 6*(x(4)/x(2))^2*B^2*R*x(2)^2*x(1) + 6*(x(4)/x(2))^2*B^2*R*ta^2*x(1) + 6*A*(x(4)/x(2))^2*B*L*x(2) + 6*A*(x(4)/x(3))^2*B*L*x(3) - 12*A*(x(4)/x(2))^2*B*L*ta - 6*(x(4)/x(2))*B*D*R*x(2)^2 - 6*(x(4)/x(3))*B*D*R*x(3)^2 - 12*(x(4)/x(2))*B*D*R*ta^2 + 6*(x(4)/x(2))*D*Ke^2*x(2)*ta + 6*(x(4)/x(2))*D*Ke^2*x(2)*x(1) - 6*(x(4)/x(2))*D*Ke^2*ta*x(1) + 6*A*(x(4)/x(2))^2*B*R*x(2)^2 + 6*A*(x(4)/x(3))^2*B*R*x(3)^2 - 6*A*(x(4)/x(2))^2*Ke^2*x(2)*ta - 6*(x(4)/x(2))^2*B^2*L*x(2)*ta + 6*(x(4)/x(2))*B*D*L*x(2) - 6*(x(4)/x(3))*B*D*L*x(3) - 12*(x(4)/x(2))*B*D*L*ta + 12*A*(x(4)/x(2))*D*R*x(2) - 12*A*(x(4)/x(3))*D*R*x(3) - 24*A*(x(4)/x(2))*D*R*ta - 6*(x(4)/x(2))*(x(4)/x(3))*B*Ke^2*x(2)*x(3)^2 + 6*(x(4)/x(2))*(x(4)/x(3))*B*Ke^2*x(3)^2*ta - 6*(x(4)/x(2))*(x(4)/x(3))*B^2*R*x(2)*x(3)^2 + 6*(x(4)/x(2))*(x(4)/x(3))*B^2*R*x(3)^2*ta - 12*(x(4)/x(2))^2*B*Ke^2*x(2)*ta*x(1) - 12*(x(4)/x(2))^2*B^2*R*x(2)*ta*x(1) + 12*(x(4)/x(2))*B*D*R*x(2)*ta + 12*(x(4)/x(2))*B*D*R*x(2)*x(1) - 12*(x(4)/x(2))*B*D*R*ta*x(1) - 6*A*(x(4)/x(2))*(x(4)/x(3))*Ke^2*x(2)*x(3) + 6*A*(x(4)/x(2))*(x(4)/x(3))*Ke^2*x(3)*ta - 6*(x(4)/x(2))*(x(4)/x(3))*B^2*L*x(2)*x(3) + 6*(x(4)/x(2))*(x(4)/x(3))*B^2*L*x(3)*ta - 12*A*(x(4)/x(2))^2*B*R*x(2)*ta - 12*A*(x(4)/x(2))*(x(4)/x(3))*B*R*x(2)*x(3) + 12*A*(x(4)/x(2))*(x(4)/x(3))*B*R*x(3)*ta)/(6*Ke^2));


options = optimoptions('ga','PlotFcn',"gaplotbestf",'Display','iter',MaxGenerations=800, ConstraintTolerance=1e-12, FunctionTolerance=1e-15);
[optxs,fval] = ga(funE, 4,[],[],[],[],[9;0.1;0.1;0.1],[9;50;50;50],nonlinfcn,options)

% options2 = optimoptions('patternsearch','PlotFcn',"psplotbestf", ConstraintTolerance=1e-12, FunctionTolerance=1e-15);
% [optxs,fval] = patternsearch(funE,[9,0.1,0.1,0.1],[],[],[],[],[9;0.1;0.1;0.1],[9;50;50;50],nonlinfcn,options2)

% options = optimoptions('particleswarm','PlotFcn',"pswplotbestf",FunctionTolerance=1e-15,MaxStallIterations=200);
% [optxs,fval] = particleswarm(funE,4,[9;0.1;0.1;0.1],[9;50;50;50],options)