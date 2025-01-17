X1 = R*D^2/Ke^2;
X2 = R*2*B*D/Ke^2 + D;
X3 = R*B^2/Ke^2 + B;
X4 = (R*A^2 + A*B*L)/Ke^2;

len = 30;

% tm = ctrl_tm;
% len = ctrl_len;

% t1 = 2;
% t3 = t1;

t1s = (0.1:0.001:20)';

% given k
tm_opt = roots([X1*(3*(1-k)/2 - 9*(1-k)^2/4 + 9*(1-k)^3/8 - 3*(1-k)^4/16) (0) X3*len^2*(5*(1-k)^2/4 - 3*(1-k)/2) (0) (12*X4*len^2)])


%tm_opt = roots([3*t1*X1 -9*t1^2*X1 ((9*t1^3*X1)-(3*X3*t1*len^2)) ((-3*t1^4*X1) + (5*t1^2*X3*len^2) + (12*X4*len^2))])
tms_opt = zeros(length(t1s), 3);
Es = zeros(length(t1s), 1);
tms = zeros(length(t1s), 1);
for i = 1:length(t1s)
    t1 = t1s(i);
    t3 = t1;
    tm_opt = roots([3*t1*X1 -9*t1^2*X1 ((9*t1^3*X1)-(3*X3*t1*len^2)) ((-3*t1^4*X1) + (5*t1^2*X3*len^2) + (12*X4*len^2))]);
%     tms_opt(i) = tm_opt(tm_opt == real(tm_opt));
    tms_opt(i,1:3) = vpa(abs(tm_opt));
    [Es(i), ind] = min(X1*tms_opt(i)  +  X2*len  +  4*len^2*X3*(tms_opt(i) - 2*t1/3 - 2*t3/3)/(2*tms_opt(i) - t1 - t3)^2  +  4*len^2*X4*(t1+t3)/((2*tms_opt(i) - t1 - t3)^2 * t1 * t3));
    tms(i) = vpa(tms_opt(i,ind));
    if(t1+t3>=tm_opt)
        break
    end
end

% figure;
% plot(tms(1:i),Es(1:i), 'b-');grid on;title('E against Tm');
% figure;
% plot(t1s(1:i),Es(1:i), 'g-');grid on;title('E against t1');
% figure;
% plot(t1s(1:i),tms(1:i), 'r-');grid on;title('Tm against t1');

% tms_opt = roots([3*t1s*X1 -9*t1s.^2*X1 ((9*t1s.^3*X1)-(3*X3*t1s*len^2)) ((-3*t1s.^4*X1) + (5*t1s.^2*X3*len^2) + (12*X4*len^2))])


% Es = zeros(size(t1s));
% for i=1:length(t1s)
%     tm = tms_opt(i);
%     t1 = t1s(i);
%     t3 = t1;
%     Es(i) = X1*tm  +  X2*len  +  4*len^2*X3*(tm - 2*t1/3 - 2*t3/3)/(2*tm - t1 - t3)^2  +  4*len^2*X4*(t1+t3)/((2*tm - t1 - t3)^2 * t1 * t3)
% end
% figure('DefaultAxesFontSize',16);
% plot(ks,Es, 'b-','LineWidth',1.5);
% title('A Graph of Energy Used Against T');
% xlabel('k');
% ylabel('Energy/J');
% grid on;
% % hold on;
% % E = X1*tm  +  X2*len  +  (2*len/(tm+tm*k))^2  * ( X3*(tm + 2*tm*k)/3 + X4*4/(tm-tm*k));
% % ks = (-6:0.02:6);
% % Es = X1*tm  +  X2*len  +  (2*len./(tm+tm*ks)).^2  .* ( X3*(tm + 2*tm*ks)/3 + X4*4./(tm-tm*ks));
% % plot(ks,Es, 'b-')
% % grid on;
% 
% k_opt = roots([X3*tm^2 -2*X3*tm^2 ((X3*tm^2)-(18*X4)) 6*X4])
% pk = poly(k_opt);
% 
% xvals=-35:0.001:39;
% 
% extrm = (8*len^2*tm^2)./(tm+tm*xvals).^3;
% % extrm = 1;
% 
% dtrms = polyval(pk,xvals);
% 
% alltrms = extrm .* dtrms;
% % alltrms = dtrms;
% % figure;
% % 
% % plot(xvals,alltrms)
% % grid on;
% 
% ta_opt = t0 + 0.5 * (tm - tm*k_opt)
% tb_opt = ta_opt + tm*k_opt
% % msk = (tb_opt > ta) & (tb_opt<tm) & (tb_opt<tc_opt);
% % k_opt = k_opt(msk);
% 
% E_opt = (X1*tm  +  X2*len  +  (2*len./(tm+tm*k_opt)).^2  .* ( X3*(tm + 2*tm*k_opt)/3 + X4*4./(tm-tm*k_opt)));
% vpa(E_opt)
% pInc =  (E-E_opt(1))*100/E