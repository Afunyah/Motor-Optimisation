X1 = R*D^2/Ke^2;
X2 = R*2*B*D/Ke^2 + D;
X3 = R*B^2/Ke^2 + B;
X4 = (R*A^2 + A*B*L)/Ke^2;

len = 30;
tm = 10;
k = 0.2;

% tm = ctrl_tm;
% len = ctrl_len;
% k = ctrl_k;

ta = 0;
tb = ta + 0.5 * (tm - k*tm);
tc = tb + k*tm;
t1 = tb;
t3 = tm - tc;

% t1 = 4.3739;
% t3 = 4.4325;

% t0 = 0;
% ta = t0 + 0.5 * (tm - k*tm);
% tb = ta + k*tm;
% t1 = ta;
% t3 = tm - tb;


% Vm = (3/(2+k)) * (len/tm);
% Vm = (2/(1+k)) * (len/tm);
% A1 = Vm/t1;
% A2 = Vm/t3;

% E expand
% E1 = X1*tm  +  X2*len  +  4*len^2*X3*(tm - 2*t1/3 - 2*t3/3)/(2*tm - t1 - t3)^2  +  4*len^2*X4*(t1+t3)/((2*tm - t1 - t3)^2 * t1 * t3)
% E2 = X1*tm  +  X2*len  +  4*len^2*X3*(tm - 2*tm/3 + 2*tm*k/3)/(2*tm - tm + tm*k)^2  +  4*len^2*X4*(tm-tm*k)/((2*tm - tm + tm*k)^2 * 0.25 * (tm-tm*k)^2)
% E3 = X1*tm  +  X2*len  +  4*len^2*X3*(tm + 2*tm*k)/(3*(tm + tm*k)^2)  +  16*len^2*X4/((tm + tm*k)^2  * (tm-tm*k))
E = X1*tm  +  X2*len  +  (2*len/(tm+tm*k))^2  * ( X3*(tm + 2*tm*k)/3 + X4*4/(tm-tm*k));

% E = (6*D^2*R*tm - 6*D^2*R*ta - 3*A1*D*Ke^2*t1^2 - 3*A2*D*Ke^2*t3^2 - 6*A1*D*Ke^2*ta^2 + 6*A^2*A1^2*R*t1 + 6*A^2*A2^2*R*t3 - 12*A^2*A1^2*R*ta + 3*A*A1^2*Ke^2*t1^2 + 3*A*A2^2*Ke^2*t3^2 - 4*A1^2*B*Ke^2*t1^3 + 2*A2^2*B*Ke^2*t3^3 + 2*A1^2*B*Ke^2*ta^3 + 3*A1^2*B^2*L*t1^2 + 3*A2^2*B^2*L*t3^2 - 4*A1^2*B^2*R*t1^3 + 2*A2^2*B^2*R*t3^3 + 2*A1^2*B^2*R*ta^3 - 12*A1^2*B*Ke^2*t1*ta^2 + 12*A1^2*B*Ke^2*t1^2*ta + 6*A1^2*B*Ke^2*t1^2*tm + 6*A1^2*B*Ke^2*ta^2*tm - 12*A1^2*B^2*R*t1*ta^2 + 12*A1^2*B^2*R*t1^2*ta + 6*A1^2*B^2*R*t1^2*tm + 6*A1^2*B^2*R*ta^2*tm + 6*A*A1^2*B*L*t1 + 6*A*A2^2*B*L*t3 - 12*A*A1^2*B*L*ta - 6*A1*B*D*R*t1^2 - 6*A2*B*D*R*t3^2 - 12*A1*B*D*R*ta^2 + 6*A1*D*Ke^2*t1*ta + 6*A1*D*Ke^2*t1*tm - 6*A1*D*Ke^2*ta*tm + 6*A*A1^2*B*R*t1^2 + 6*A*A2^2*B*R*t3^2 - 6*A*A1^2*Ke^2*t1*ta - 6*A1^2*B^2*L*t1*ta + 6*A1*B*D*L*t1 - 6*A2*B*D*L*t3 - 12*A1*B*D*L*ta + 12*A*A1*D*R*t1 - 12*A*A2*D*R*t3 - 24*A*A1*D*R*ta - 6*A1*A2*B*Ke^2*t1*t3^2 + 6*A1*A2*B*Ke^2*t3^2*ta - 6*A1*A2*B^2*R*t1*t3^2 + 6*A1*A2*B^2*R*t3^2*ta - 12*A1^2*B*Ke^2*t1*ta*tm - 12*A1^2*B^2*R*t1*ta*tm + 12*A1*B*D*R*t1*ta + 12*A1*B*D*R*t1*tm - 12*A1*B*D*R*ta*tm - 6*A*A1*A2*Ke^2*t1*t3 + 6*A*A1*A2*Ke^2*t3*ta - 6*A1*A2*B^2*L*t1*t3 + 6*A1*A2*B^2*L*t3*ta - 12*A*A1^2*B*R*t1*ta - 12*A*A1*A2*B*R*t1*t3 + 12*A*A1*A2*B*R*t3*ta)/(6*Ke^2);
% E = -(20*A^2*L*Vm^2*t1^4 - 20*A^2*R*Vm^2*t1^4 + 15*D^2*R*t1^4*t3*ta - 15*D^2*R*t1^4*t3*tm + 20*A^2*L*Vm^2*t1^3*t3 + 20*A^2*L*Vm^2*t3*ta^3 + 7*B*Ke^2*Vm^2*t1^5*t3 + 6*B*Ke^2*Vm^2*t3*ta^5 + 5*D*Ke^2*Vm*t1^4*t3^2 - 20*A^2*R*Vm^2*t1^3*t3 + 40*A^2*R*Vm^2*t3*ta^3 + 7*B^2*R*Vm^2*t1^5*t3 + 6*B^2*R*Vm^2*t3*ta^5 + 30*A*D*L*Vm*t1^4 + 7*B*Ke^2*Vm^2*t1^4*t3^2 + 7*B^2*R*Vm^2*t1^4*t3^2 + 5*D*Ke^2*Vm*t1^5*t3 + 45*A*Ke^2*Vm^2*t1^2*t3*ta^2 + 10*B*Ke^2*Vm^2*t1^2*t3*ta^3 + 45*B^2*L*Vm^2*t1^2*t3*ta^2 + 10*B^2*R*Vm^2*t1^2*t3*ta^3 + 60*A*B*L*Vm^2*t3*ta^3 + 10*B*D*R*Vm*t1^4*t3^2 - 15*D*Ke^2*Vm*t1^4*t3*tm - 30*A*Ke^2*Vm^2*t1*t3*ta^3 - 30*A^2*L*Vm^2*t1*t3*ta^2 - 30*A^2*L*Vm^2*t1^2*t3*ta - 15*B*Ke^2*Vm^2*t1*t3*ta^4 - 15*B*Ke^2*Vm^2*t1^4*t3*tm - 30*B^2*L*Vm^2*t1*t3*ta^3 - 10*D*Ke^2*Vm*t1^2*t3*ta^3 + 15*D*Ke^2*Vm*t1^3*t3*ta^2 - 60*A^2*R*Vm^2*t1*t3*ta^2 + 60*A^2*R*Vm^2*t1^2*t3*ta - 15*B^2*R*Vm^2*t1*t3*ta^4 - 15*B^2*R*Vm^2*t1^4*t3*tm + 30*A*D*L*Vm*t1^3*t3 + 10*B*D*R*Vm*t1^5*t3 + 90*A*B*R*Vm^2*t1^2*t3*ta^2 - 60*A*D*L*Vm*t1^2*t3*ta + 30*B*D*L*Vm*t1^3*t3*ta + 60*A*D*R*Vm*t1^3*t3*ta - 30*B*D*R*Vm*t1^4*t3*tm - 90*A*B*L*Vm^2*t1*t3*ta^2 + 30*A*B*L*Vm^2*t1^2*t3*ta - 60*A*B*R*Vm^2*t1*t3*ta^3 - 20*B*D*R*Vm*t1^2*t3*ta^3 + 30*B*D*R*Vm*t1^3*t3*ta^2)/(15*Ke^2*t1^4*t3);
vpa(E)


ks = (0:0.0001:0.9);
Es = zeros(size(ks));
for i=1:length(ks)
    Es(i) = X1*tm  +  X2*len  +  (2*len/(tm+tm*ks(i)))^2  * ( X3*(tm + 2*tm*ks(i))/3 + X4*4/(tm-tm*ks(i)));
k = ks(i);
ta = 0;
tb = ta + 0.5 * (tm - k*tm);
tc = tb + k*tm;
t1 = tb;
t3 = tm - tc;
% Vm = (3/(2+k)) * (len/tm);
%     Es(i) = -(20*A^2*L*Vm^2*t1^4 - 20*A^2*R*Vm^2*t1^4 + 15*D^2*R*t1^4*t3*ta - 15*D^2*R*t1^4*t3*tm + 20*A^2*L*Vm^2*t1^3*t3 + 20*A^2*L*Vm^2*t3*ta^3 + 7*B*Ke^2*Vm^2*t1^5*t3 + 6*B*Ke^2*Vm^2*t3*ta^5 + 5*D*Ke^2*Vm*t1^4*t3^2 - 20*A^2*R*Vm^2*t1^3*t3 + 40*A^2*R*Vm^2*t3*ta^3 + 7*B^2*R*Vm^2*t1^5*t3 + 6*B^2*R*Vm^2*t3*ta^5 + 30*A*D*L*Vm*t1^4 + 7*B*Ke^2*Vm^2*t1^4*t3^2 + 7*B^2*R*Vm^2*t1^4*t3^2 + 5*D*Ke^2*Vm*t1^5*t3 + 45*A*Ke^2*Vm^2*t1^2*t3*ta^2 + 10*B*Ke^2*Vm^2*t1^2*t3*ta^3 + 45*B^2*L*Vm^2*t1^2*t3*ta^2 + 10*B^2*R*Vm^2*t1^2*t3*ta^3 + 60*A*B*L*Vm^2*t3*ta^3 + 10*B*D*R*Vm*t1^4*t3^2 - 15*D*Ke^2*Vm*t1^4*t3*tm - 30*A*Ke^2*Vm^2*t1*t3*ta^3 - 30*A^2*L*Vm^2*t1*t3*ta^2 - 30*A^2*L*Vm^2*t1^2*t3*ta - 15*B*Ke^2*Vm^2*t1*t3*ta^4 - 15*B*Ke^2*Vm^2*t1^4*t3*tm - 30*B^2*L*Vm^2*t1*t3*ta^3 - 10*D*Ke^2*Vm*t1^2*t3*ta^3 + 15*D*Ke^2*Vm*t1^3*t3*ta^2 - 60*A^2*R*Vm^2*t1*t3*ta^2 + 60*A^2*R*Vm^2*t1^2*t3*ta - 15*B^2*R*Vm^2*t1*t3*ta^4 - 15*B^2*R*Vm^2*t1^4*t3*tm + 30*A*D*L*Vm*t1^3*t3 + 10*B*D*R*Vm*t1^5*t3 + 90*A*B*R*Vm^2*t1^2*t3*ta^2 - 60*A*D*L*Vm*t1^2*t3*ta + 30*B*D*L*Vm*t1^3*t3*ta + 60*A*D*R*Vm*t1^3*t3*ta - 30*B*D*R*Vm*t1^4*t3*tm - 90*A*B*L*Vm^2*t1*t3*ta^2 + 30*A*B*L*Vm^2*t1^2*t3*ta - 60*A*B*R*Vm^2*t1*t3*ta^3 - 20*B*D*R*Vm*t1^2*t3*ta^3 + 30*B*D*R*Vm*t1^3*t3*ta^2)/(15*Ke^2*t1^4*t3);

end
figure('DefaultAxesFontSize',16);
plot(ks,Es, 'b-','LineWidth',1.5);
title('A Graph of Energy Used Against k');
xlabel('k');
ylabel('Energy/J');
grid on;
% hold on;
% E = X1*tm  +  X2*len  +  (2*len/(tm+tm*k))^2  * ( X3*(tm + 2*tm*k)/3 + X4*4/(tm-tm*k));
% ks = (-6:0.02:6);
% Es = X1*tm  +  X2*len  +  (2*len./(tm+tm*ks)).^2  .* ( X3*(tm + 2*tm*ks)/3 + X4*4./(tm-tm*ks));
% plot(ks,Es, 'b-')
% grid on;

k_opt = roots([X3*tm^2 -2*X3*tm^2 ((X3*tm^2)-(18*X4)) 6*X4]);
pk = poly(k_opt);

xvals=-35:0.001:39;

extrm = (8*len^2*tm^2)./(tm+tm*xvals).^3;
% extrm = 1;

dtrms = polyval(pk,xvals);

alltrms = extrm .* dtrms;
% alltrms = dtrms;
% figure;
% 
% plot(xvals,alltrms)
% grid on;

% ta_opt = t0 + 0.5 * (tm - tm*k_opt)
% tb_opt = ta_opt + tm*k_opt
% 
% 
% E_opt = (X1*tm  +  X2*len  +  (2*len./(tm+tm*k_opt)).^2  .* ( X3*(tm + 2*tm*k_opt)/3 + X4*4./(tm-tm*k_opt)));
% vpa(E_opt)
% pInc =  (E-E_opt(1))*100/E