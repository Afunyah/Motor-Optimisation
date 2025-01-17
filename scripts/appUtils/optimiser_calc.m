USE_ACC_K = 0;
USE_ACC_TM = 0;
USE_ACC = 0;

USE_VEL_K = 0;
USE_VEL_TM = 0;
USE_VEL = 0;

USE_ACC_VEL = 0;

USE_K = 0;
USE_TM = 0;

USE_K_TM = 0;

NONE = 0;

switch(opt_params)
    case 'NONE'; NONE = 1;
    case 'tm,t1,t3,Vm - given acc,dec,k'; USE_ACC_K = 1;
    case 'k,t1,t3,Vm - given acc,dec,tm'; USE_ACC_TM = 1;
    case 'ALL - given acc,dec'; USE_ACC = 1;
    case 'tm,t1,t3,acc,dec - given Vm,k'; USE_VEL_K = 1;
    case 'k,t1,t3,acc,dec - given Vm,tm'; USE_VEL_TM = 1;
    case 'ALL - given Vm'; USE_VEL = 1;
    case 'tm,k,t1,t3 - given Vm,acc,dec'; USE_ACC_VEL = 1;
    case 'ALL - given k'; USE_K = 1;
    case 'ALL - given tm'; USE_TM = 1;
    case 't1,t3,Vm,acc,dec - given k,tm'; USE_K_TM = 1;
    otherwise; error('controller/optimiser calculator error');
end

if(USE_ACC_K+USE_ACC_TM+USE_ACC+USE_VEL_K+USE_VEL_TM+USE_VEL+USE_ACC_VEL+USE_K+USE_TM+USE_K_TM+NONE ~= 1)
    error('controller/optimiser calculator error');
end

mp = 0;
if(strcmpi(MP,'trapezoidal') == 1)
    mp = 1;
    funE = @(x) ((6*D^2*R*x(1) - 6*D^2*R*ta - 3*(x(4)/x(2))*D*Ke^2*x(2)^2 - 3*(x(4)/x(3))*D*Ke^2*x(3)^2 - 6*(x(4)/x(2))*D*Ke^2*ta^2 + 6*A^2*(x(4)/x(2))^2*R*x(2) + 6*A^2*(x(4)/x(3))^2*R*x(3) - 12*A^2*(x(4)/x(2))^2*R*ta + 3*A*(x(4)/x(2))^2*Ke^2*x(2)^2 + 3*A*(x(4)/x(3))^2*Ke^2*x(3)^2 - 4*(x(4)/x(2))^2*B*Ke^2*x(2)^3 + 2*(x(4)/x(3))^2*B*Ke^2*x(3)^3 + 2*(x(4)/x(2))^2*B*Ke^2*ta^3 + 3*(x(4)/x(2))^2*B^2*L*x(2)^2 + 3*(x(4)/x(3))^2*B^2*L*x(3)^2 - 4*(x(4)/x(2))^2*B^2*R*x(2)^3 + 2*(x(4)/x(3))^2*B^2*R*x(3)^3 + 2*(x(4)/x(2))^2*B^2*R*ta^3 - 12*(x(4)/x(2))^2*B*Ke^2*x(2)*ta^2 + 12*(x(4)/x(2))^2*B*Ke^2*x(2)^2*ta + 6*(x(4)/x(2))^2*B*Ke^2*x(2)^2*x(1) + 6*(x(4)/x(2))^2*B*Ke^2*ta^2*x(1) - 12*(x(4)/x(2))^2*B^2*R*x(2)*ta^2 + 12*(x(4)/x(2))^2*B^2*R*x(2)^2*ta + 6*(x(4)/x(2))^2*B^2*R*x(2)^2*x(1) + 6*(x(4)/x(2))^2*B^2*R*ta^2*x(1) + 6*A*(x(4)/x(2))^2*B*L*x(2) + 6*A*(x(4)/x(3))^2*B*L*x(3) - 12*A*(x(4)/x(2))^2*B*L*ta - 6*(x(4)/x(2))*B*D*R*x(2)^2 - 6*(x(4)/x(3))*B*D*R*x(3)^2 - 12*(x(4)/x(2))*B*D*R*ta^2 + 6*(x(4)/x(2))*D*Ke^2*x(2)*ta + 6*(x(4)/x(2))*D*Ke^2*x(2)*x(1) - 6*(x(4)/x(2))*D*Ke^2*ta*x(1) + 6*A*(x(4)/x(2))^2*B*R*x(2)^2 + 6*A*(x(4)/x(3))^2*B*R*x(3)^2 - 6*A*(x(4)/x(2))^2*Ke^2*x(2)*ta - 6*(x(4)/x(2))^2*B^2*L*x(2)*ta + 6*(x(4)/x(2))*B*D*L*x(2) - 6*(x(4)/x(3))*B*D*L*x(3) - 12*(x(4)/x(2))*B*D*L*ta + 12*A*(x(4)/x(2))*D*R*x(2) - 12*A*(x(4)/x(3))*D*R*x(3) - 24*A*(x(4)/x(2))*D*R*ta - 6*(x(4)/x(2))*(x(4)/x(3))*B*Ke^2*x(2)*x(3)^2 + 6*(x(4)/x(2))*(x(4)/x(3))*B*Ke^2*x(3)^2*ta - 6*(x(4)/x(2))*(x(4)/x(3))*B^2*R*x(2)*x(3)^2 + 6*(x(4)/x(2))*(x(4)/x(3))*B^2*R*x(3)^2*ta - 12*(x(4)/x(2))^2*B*Ke^2*x(2)*ta*x(1) - 12*(x(4)/x(2))^2*B^2*R*x(2)*ta*x(1) + 12*(x(4)/x(2))*B*D*R*x(2)*ta + 12*(x(4)/x(2))*B*D*R*x(2)*x(1) - 12*(x(4)/x(2))*B*D*R*ta*x(1) - 6*A*(x(4)/x(2))*(x(4)/x(3))*Ke^2*x(2)*x(3) + 6*A*(x(4)/x(2))*(x(4)/x(3))*Ke^2*x(3)*ta - 6*(x(4)/x(2))*(x(4)/x(3))*B^2*L*x(2)*x(3) + 6*(x(4)/x(2))*(x(4)/x(3))*B^2*L*x(3)*ta - 12*A*(x(4)/x(2))^2*B*R*x(2)*ta - 12*A*(x(4)/x(2))*(x(4)/x(3))*B*R*x(2)*x(3) + 12*A*(x(4)/x(2))*(x(4)/x(3))*B*R*x(3)*ta)/(6*Ke^2));
elseif(strcmpi(MP,'cosine') == 1)
    mp = 2;
    funE = @(x) (-(2*A^2*L*x(4)^2*x(2)^2*pi^3 - 2*A^2*R*x(4)^2*x(2)^2*pi^3 + 8*pi*D*Ke^2*x(4)*x(2)^3*x(3) - 8*D*Ke^2*x(4)*x(2)^3*x(3)*sin((pi*ta)/x(2)) + 16*pi*D^2*R*x(2)^2*x(3)*ta - 16*pi*D^2*R*x(2)^2*x(3)*x(1) + 8*D*Ke^2*x(4)*x(2)^3*x(3)*sin((pi*(x(2) - ta))/x(2)) + 8*pi*A*Ke^2*x(4)^2*x(2)^2*x(3) + 2*A^2*L*x(4)^2*x(2)*x(3)*pi^3 - 4*A^2*L*x(4)^2*x(3)*ta*pi^3 + 10*pi*B*Ke^2*x(4)^2*x(2)^3*x(3) + 8*pi*B^2*L*x(4)^2*x(2)^2*x(3) + 8*pi*D*Ke^2*x(4)*x(2)^2*x(3)^2 - 2*A^2*R*x(4)^2*x(2)*x(3)*pi^3 + 4*A^2*R*x(4)^2*x(3)*ta*pi^3 + 10*pi*B^2*R*x(4)^2*x(2)^3*x(3) - 8*B*Ke^2*x(4)^2*x(2)^3*x(3)*sin((pi*ta)/x(2)) + B*Ke^2*x(4)^2*x(2)^3*x(3)*sin((2*pi*ta)/x(2)) - 8*B^2*R*x(4)^2*x(2)^3*x(3)*sin((pi*ta)/x(2)) + B^2*R*x(4)^2*x(2)^3*x(3)*sin((2*pi*ta)/x(2)) + 8*B*Ke^2*x(4)^2*x(2)^3*x(3)*sin((pi*(x(2) - ta))/x(2)) - B*Ke^2*x(4)^2*x(2)^3*x(3)*sin((2*pi*(x(2) - ta))/x(2)) + 8*B^2*R*x(4)^2*x(2)^3*x(3)*sin((pi*(x(2) - ta))/x(2)) - B^2*R*x(4)^2*x(2)^3*x(3)*sin((2*pi*(x(2) - ta))/x(2)) + 10*pi*B*Ke^2*x(4)^2*x(2)^2*x(3)^2 + 10*pi*B^2*R*x(4)^2*x(2)^2*x(3)^2 + 4*A^2*L*x(4)^2*x(2)*x(3)*pi^2*sin((pi*ta)/x(2)) - A^2*L*x(4)^2*x(2)*x(3)*pi^2*sin((2*pi*ta)/x(2)) - A^2*R*x(4)^2*x(2)*x(3)*pi^2*sin((2*pi*ta)/x(2)) - 4*pi*B*Ke^2*x(4)^2*x(2)^2*x(3)*ta - 16*pi*B*Ke^2*x(4)^2*x(2)^2*x(3)*x(1) - 4*pi*B^2*R*x(4)^2*x(2)^2*x(3)*ta - 16*pi*B^2*R*x(4)^2*x(2)^2*x(3)*x(1) + 16*pi*B*D*L*x(4)*x(2)^2*x(3) + 32*pi*A*D*R*x(4)*x(2)^2*x(3) + 16*pi*B*D*R*x(4)*x(2)^3*x(3) + 4*A*Ke^2*x(4)^2*x(2)^2*x(3)*pi*cos((pi*(x(2) - ta))/x(2)) - A*Ke^2*x(4)^2*x(2)^2*x(3)*pi*cos((2*pi*(x(2) - ta))/x(2)) + 4*B^2*L*x(4)^2*x(2)^2*x(3)*pi*cos((pi*(x(2) - ta))/x(2)) - B^2*L*x(4)^2*x(2)^2*x(3)*pi*cos((2*pi*(x(2) - ta))/x(2)) - 4*A^2*L*x(4)^2*x(2)*x(3)*pi^2*sin((pi*(x(2) - ta))/x(2)) + A^2*L*x(4)^2*x(2)*x(3)*pi^2*sin((2*pi*(x(2) - ta))/x(2)) + A^2*R*x(4)^2*x(2)*x(3)*pi^2*sin((2*pi*(x(2) - ta))/x(2)) - 16*B*D*R*x(4)*x(2)^3*x(3)*sin((pi*ta)/x(2)) + 16*B*D*R*x(4)*x(2)^3*x(3)*sin((pi*(x(2) - ta))/x(2)) + 16*pi*A*B*R*x(4)^2*x(2)^2*x(3) + 16*pi*B*D*R*x(4)*x(2)^2*x(3)^2 - 16*pi*D*Ke^2*x(4)*x(2)^2*x(3)*x(1) - 4*A*Ke^2*x(4)^2*x(2)^2*x(3)*pi*cos((pi*ta)/x(2)) + A*Ke^2*x(4)^2*x(2)^2*x(3)*pi*cos((2*pi*ta)/x(2)) - 4*B^2*L*x(4)^2*x(2)^2*x(3)*pi*cos((pi*ta)/x(2)) + B^2*L*x(4)^2*x(2)^2*x(3)*pi*cos((2*pi*ta)/x(2)) + 8*A*B*R*x(4)^2*x(2)^2*x(3)*pi*cos((pi*(x(2) - ta))/x(2)) - 2*A*B*R*x(4)^2*x(2)^2*x(3)*pi*cos((2*pi*(x(2) - ta))/x(2)) - 4*A*B*L*x(4)^2*x(2)*x(3)*pi^2*sin((pi*(x(2) - ta))/x(2)) + 2*A*B*L*x(4)^2*x(2)*x(3)*pi^2*sin((2*pi*(x(2) - ta))/x(2)) - 8*B*D*L*x(4)*x(2)^2*x(3)*pi*cos((pi*ta)/x(2)) - 16*A*D*R*x(4)*x(2)^2*x(3)*pi*cos((pi*ta)/x(2)) + 8*A*D*L*x(4)*x(2)*x(3)*pi^2*sin((pi*ta)/x(2)) - 32*pi*B*D*R*x(4)*x(2)^2*x(3)*x(1) + 8*B*D*L*x(4)*x(2)^2*x(3)*pi*cos((pi*(x(2) - ta))/x(2)) + 16*A*D*R*x(4)*x(2)^2*x(3)*pi*cos((pi*(x(2) - ta))/x(2)) - 8*A*D*L*x(4)*x(2)*x(3)*pi^2*sin((pi*(x(2) - ta))/x(2)) - 8*A*B*R*x(4)^2*x(2)^2*x(3)*pi*cos((pi*ta)/x(2)) + 2*A*B*R*x(4)^2*x(2)^2*x(3)*pi*cos((2*pi*ta)/x(2)) + 4*A*B*L*x(4)^2*x(2)*x(3)*pi^2*sin((pi*ta)/x(2)) - 2*A*B*L*x(4)^2*x(2)*x(3)*pi^2*sin((2*pi*ta)/x(2)))/(16*Ke^2*x(2)^2*x(3)*pi));
elseif(strcmpi(MP,'cubic') == 1)
    mp = 3;
    funE = @(x) (-(20*A^2*L*x(4)^2*x(2)^4 - 20*A^2*R*x(4)^2*x(2)^4 + 15*D^2*R*x(2)^4*x(3)*ta - 15*D^2*R*x(2)^4*x(3)*x(1) + 20*A^2*L*x(4)^2*x(2)^3*x(3) + 20*A^2*L*x(4)^2*x(3)*ta^3 + 7*B*Ke^2*x(4)^2*x(2)^5*x(3) + 6*B*Ke^2*x(4)^2*x(3)*ta^5 + 5*D*Ke^2*x(4)*x(2)^4*x(3)^2 - 20*A^2*R*x(4)^2*x(2)^3*x(3) + 40*A^2*R*x(4)^2*x(3)*ta^3 + 7*B^2*R*x(4)^2*x(2)^5*x(3) + 6*B^2*R*x(4)^2*x(3)*ta^5 + 30*A*D*L*x(4)*x(2)^4 + 7*B*Ke^2*x(4)^2*x(2)^4*x(3)^2 + 7*B^2*R*x(4)^2*x(2)^4*x(3)^2 + 5*D*Ke^2*x(4)*x(2)^5*x(3) + 45*A*Ke^2*x(4)^2*x(2)^2*x(3)*ta^2 + 10*B*Ke^2*x(4)^2*x(2)^2*x(3)*ta^3 + 45*B^2*L*x(4)^2*x(2)^2*x(3)*ta^2 + 10*B^2*R*x(4)^2*x(2)^2*x(3)*ta^3 + 60*A*B*L*x(4)^2*x(3)*ta^3 + 10*B*D*R*x(4)*x(2)^4*x(3)^2 - 15*D*Ke^2*x(4)*x(2)^4*x(3)*x(1) - 30*A*Ke^2*x(4)^2*x(2)*x(3)*ta^3 - 30*A^2*L*x(4)^2*x(2)*x(3)*ta^2 - 30*A^2*L*x(4)^2*x(2)^2*x(3)*ta - 15*B*Ke^2*x(4)^2*x(2)*x(3)*ta^4 - 15*B*Ke^2*x(4)^2*x(2)^4*x(3)*x(1) - 30*B^2*L*x(4)^2*x(2)*x(3)*ta^3 - 10*D*Ke^2*x(4)*x(2)^2*x(3)*ta^3 + 15*D*Ke^2*x(4)*x(2)^3*x(3)*ta^2 - 60*A^2*R*x(4)^2*x(2)*x(3)*ta^2 + 60*A^2*R*x(4)^2*x(2)^2*x(3)*ta - 15*B^2*R*x(4)^2*x(2)*x(3)*ta^4 - 15*B^2*R*x(4)^2*x(2)^4*x(3)*x(1) + 30*A*D*L*x(4)*x(2)^3*x(3) + 10*B*D*R*x(4)*x(2)^5*x(3) + 90*A*B*R*x(4)^2*x(2)^2*x(3)*ta^2 - 60*A*D*L*x(4)*x(2)^2*x(3)*ta + 30*B*D*L*x(4)*x(2)^3*x(3)*ta + 60*A*D*R*x(4)*x(2)^3*x(3)*ta - 30*B*D*R*x(4)*x(2)^4*x(3)*x(1) - 90*A*B*L*x(4)^2*x(2)*x(3)*ta^2 + 30*A*B*L*x(4)^2*x(2)^2*x(3)*ta - 60*A*B*R*x(4)^2*x(2)*x(3)*ta^3 - 20*B*D*R*x(4)*x(2)^2*x(3)*ta^3 + 30*B*D*R*x(4)*x(2)^3*x(3)*ta^2)/(15*Ke^2*x(2)^4*x(3)));
end


if(USE_ACC_K)
    c = @(x)[];
    switch(mp)
        case 1; ceq = @(x)[2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; (x(1)-x(2)-x(3))/x(1) - k; x(2)+x(3)+k*x(1)-x(1); x(4)/max_acc-x(2); x(4)/max_dec-x(3);];
        case 2; ceq = @(x)[2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; (x(1)-x(2)-x(3))/x(1) - k; x(4)*pi/(2*x(2))-max_acc; x(4)*pi/(2*x(3))-max_dec; (2/(1+k))*(len/x(1))-x(4)];
        case 3; ceq = @(x)[3*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-2; (x(1)-x(2)-x(3))/x(1) - k; 2*x(4)/x(2)-max_acc; 2*x(4)/x(3)-max_dec; ];
    end
    nonlinfcn = @(x)deal(c(x),ceq(x));
    lb = [0.1;0.1;0.1;0.1];
    ub = [50;50;50;50];
elseif(USE_ACC_TM)
    
    switch(mp)
        case 1
            c = @(x) [x(2)+x(3)-x(1); x(4)/max_acc-x(2); x(4)/max_dec-x(3);]; 
            ceq = @(x)[2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1;];
        case 2 
            c = @(x) [x(2)+x(3)-x(1); x(4)*pi/(2*x(2))-max_acc; x(4)*pi/(2*x(3))-max_dec;]; 
            ceq = @(x)[2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1;  x(4)*(x(1)/2-x(2)/2)-len/2;];
        case 3
            c = @(x) [x(2)+x(3)-x(1); x(4)*2/x(2)-max_acc; x(4)*2/x(3)-max_dec;]; 
            ceq = @(x)[3*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-2;  2*x(4)*(x(2)+x(3))/3+x(4)*(x(1)-x(2)-x(3))-len;];
    end
    nonlinfcn = @(x)deal(c(x),ceq(x));
    lb = [tm;0.1;0.1;0.1];
    ub = [tm;50;50;50];
elseif(USE_ACC)
    c = @(x) x(2)+x(3)-x(1);
    switch(mp)
        case 1; ceq = @(x)[2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; x(4)/max_acc-x(2); x(4)/max_dec-x(3);];
        case 2; ceq = @(x)[2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; x(4)*pi/(2*x(2))-max_acc; x(4)*pi/(2*x(3))-max_dec;];
        case 3; ceq = @(x)[3*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-2; x(4)*2/x(2)-max_acc; x(4)*2/x(3)-max_dec; ];
    end
    nonlinfcn = @(x)deal(c(x),ceq(x));
    lb = [0.1;0.1;0.1;0.1];
    ub = [50;50;50;50];
elseif(USE_VEL_K)
    c = @(x) x(2)+x(3)-x(1);
    switch(mp)
        case 1; ceq = @(x)[2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; (x(1)-x(2)-x(3))/x(1) - k; (2/(1+k))*(len/x(1))-x(4)];
        case 2; ceq = @(x)[2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; (x(1)-x(2)-x(3))/x(1) - k;];
        case 3; ceq = @(x)[3*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-2; (x(1)-x(2)-x(3))/x(1) - k;];
    end
    nonlinfcn = @(x)deal(c(x),ceq(x));
    lb = [0.1;0.1;0.1;max_vel];
    ub = [50;50;50;max_vel];
elseif(USE_VEL_TM)
    if(max_vel <= len/tm)
        error('MAX_VEL too low');
    end
    c = @(x) x(2)+x(3)-x(1);
    switch(mp)
        case 1; ceq = @(x)[2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1;  len-x(4)*((2*x(1)-x(2)-x(3))/2)];
        case 2; ceq = @(x)[2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; x(4)*(x(1)-x(2)/2-x(3))+x(4)*x(3)/2-len;];
        case 3; ceq = @(x)[3*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-2; 2*x(4)*(x(2)+x(3))/3+x(4)*(x(1)-x(2)-x(3))-len;];
    end
    nonlinfcn = @(x)deal(c(x),ceq(x));
    lb = [0.1;0.1;0.1;0.1];
    ub = [tm;50;50;max_vel];
elseif(USE_VEL)
    c = @(x) x(2)+x(3)-x(1);
    switch(mp)
        case 1; ceq = @(x)2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1;
        case 2; ceq = @(x)2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1;
        case 3; ceq = @(x)3*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-2;
    end
    nonlinfcn = @(x)deal(c(x),ceq(x));
    lb = [0.1;0.1;0.1;max_vel];
    ub = [50;50;50;max_vel];
elseif(USE_ACC_VEL)
    c = @(x) [x(2)+x(3)-x(1); x(4)*2/x(2)-max_acc; x(4)*2/x(3)-max_dec;];

    switch(mp)
        case 1; ceq = @(x) 2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1;
        case 2; ceq = @(x) 2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1;
        case 3; ceq = @(x) 3*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-2;
    end
    nonlinfcn = @(x)deal(c(x),ceq(x));
    lb = [0.1;0.1;0.1;0.1];
    ub = [50;50;50;max_vel];
elseif(USE_K)
    c = @(x) x(2)+x(3)-x(1);
    switch(mp)
        case 1; ceq = @(x) [2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; (x(1)-x(2)-x(3))/x(1) - k; x(2)+x(3)+k*x(1)-x(1); (2/(1+k))*(len/x(1))-x(4);];
        case 2; ceq = @(x) [2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; (x(1)-x(2)-x(3))/x(1) - k; x(2)+x(3)+k*x(1)-x(1); (2/(1+k))*(len/x(1))-x(4);];
        case 3; ceq = @(x) [3*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-2; (x(1)-x(2)-x(3))/x(1) - k; x(2)+x(3)+k*x(1)-x(1); (3/(2+k))*(len/x(1))-x(4);];
    end
    nonlinfcn = @(x)deal(c(x),ceq(x));
    lb = [0.1;0.1;0.1;0.1];
    ub = [50;50;50;50];
elseif(USE_TM)
    c = @(x) x(2)+x(3)-x(1);
    ceq = @(x)2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1;
    switch(mp)
        case 1; ceq = @(x) 2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1;
        case 2; ceq = @(x) 2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1;
        case 3; ceq = @(x) 3*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-2;
    end
    nonlinfcn = @(x)deal(c(x),ceq(x));
    lb = [tm;0.1;0.1;0.1];
    ub = [tm;50;50;50];
elseif(USE_K_TM)
    c = @(x) x(2)+x(3)-x(1);
    ceq = @(x) [2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; (x(1)-x(2)-x(3))/x(1) - k; x(2)+x(3)+k*x(1)-x(1);];
    switch(mp)
        case 1; ceq = @(x) [2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; (x(1)-x(2)-x(3))/x(1) - k; x(2)+x(3)+k*x(1)-x(1);];
        case 2; ceq = @(x) [2*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-1; (x(1)-x(2)-x(3))/x(1) - k; x(2)+x(3)+k*x(1)-x(1);];
        case 3; ceq = @(x) [3*len/(x(1)*x(4))-(x(1)-x(2)-x(3))/x(1)-2; (x(1)-x(2)-x(3))/x(1) - k; x(2)+x(3)+k*x(1)-x(1);];
    end
    nonlinfcn = @(x)deal(c(x),ceq(x));
    lb = [tm;0.1;0.1;0.1];
    ub = [tm;50;50;50];
end

if(NONE)
    switch(lower(no_opt_params))
        case 't1 and t3'
        case 'max acceleration and deceleration'
            switch(mp)
                case 1
                    t1 = max_vel/max_acc;
                    t3 = max_vel/max_dec;
                case 2
                    t1 = max_vel*pi/(2*max_acc);
                    t3 = max_vel*pi/(2*max_dec);
                case 3
                    t1 = max_vel*2/(max_acc);
                    t3 = max_vel*2/(max_dec);
            end
    end
    E_total = funE([tm,t1,t3,max_vel]);
    k = (tm-t1-t3)/tm;
    cmes1 = '';cmes2 = '';
else
    hybridopts = optimoptions('fmincon','Display','none','Algorithm','interior-point','OutputFcn',@optimplotfval1,EnableFeasibilityMode=true,ConstraintTolerance=1e-7, OptimalityTolerance=1e-8);
    options = optimoptions('ga','OutputFcn',@(options,state,flag)gaplotbestf1(options,state,flag),SelectionFcn="selectionstochunif",MutationFcn="mutationadaptfeasible",MaxStallGenerations=35,MaxGenerations=500, ConstraintTolerance=1e-6, FunctionTolerance=1e-14);
    
    [ga_OPTXS,ga_FVAL,ga_EXFLAG,OUTPUT] = ga(funE, 4,[],[],[],[],lb,ub, nonlinfcn,options);
    cmes1 = OUTPUT.message;
    use_ga_res = 0;

    [OPTXS,FVAL,fmi_EXFLAG,OUTPUT] = fmincon(funE,ga_OPTXS,[],[],[],[],lb,ub,nonlinfcn,hybridopts);
    cmes2 = OUTPUT.message;
    
    if(fmi_EXFLAG == -2 && ga_EXFLAG ~= -2)
        OPTXS = ga_OPTXS;
        FVAL = ga_FVAL;
        use_ga_res = 1;
    end

    tm = OPTXS(1); t1 = OPTXS(2); t3 = OPTXS(3); max_vel = OPTXS(4);
    k = (tm-t1-t3)/tm;
    E_total = FVAL;
end

calc_opt_res = sprintf('tm = %0.2f s \nt1 = %0.2f s \nt3 = %0.2f s \nmax_vel = %0.2f rad/s \nk = %0.2f\nE = %0.3f J',tm,t1,t3,max_vel,k,E_total);
if(use_ga_res)
    calc_opt_message = cat(2,cmes1,cmes2,calc_opt_res,"\nUSING GA RESULTS\n"); %#ok<*NASGU> 
end
calc_opt_message = cat(2,cmes1,cmes2,calc_opt_res);

run("controller_calc.m");
