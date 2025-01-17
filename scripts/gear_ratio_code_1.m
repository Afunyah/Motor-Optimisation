
A1 = J_m + J_s/n_m + J_g/n_m;
A2 = (2*J_p + J_l + J_b + J_s)/(n_m*n_g);

B1 = b_m + b/n_m ;
B2 = b/(n_m*n_g);

D1 = (mu_s*cos(angl) + sin(angl)) * r_p*(m_l+m_b)*g/(n_m*n_g);

X11 = (D1^2*R)/(Ke^2);

X21 = (2*B2*D1*R)/(Ke^2);
X22 = (2*B1*D1*R)/(Ke^2) + D1;

X31 = (B2^2*R)/(Ke^2);
X32 = B2 + (2*B1*B2*R)/(Ke^2);
X33 = (B1^2*R)/Ke^2 + B1;

X41 = (A2^2*R)/(Ke^2) + (A2*B2*L)/(Ke^2);
X42 = (A1*B2*L + A2*B1*L + 2*A1*A2*R)/Ke^2;
X43 = (A1^2 + A1*B1*L)/Ke^2;

myfun = @(N,X11,X21,X22,X31,X32,X41,X42,tm,len,t1,t3) (-2*X11*tm)/N^3 - (3*X21*len)/N^4 - (X22*len)/N^2 + (-(4*X31)/(N^5) - (2*X32)/(N^3))*(4*len^2*(tm - 2*t1/3 - 2*t3/3)/(2*tm - t1 - t3)^2) + (-(4*X41)/(N^5) - (2*X42)/(N^3))*(4*len^2*(t1+t3)/((2*tm - t1 - t3)^2 * t1 * t3));
dEdN = @(N) myfun(N,X11,X21,X22,X31,X32,X41,X42,tm,len,t1,t3);
% options = optimset('PlotFcns',{@optimplotx,@optimplotfval});

p = fzero(dEdN,1)

N = 1;
E = (X11*tm)/N^2 + (X21*len)/N^3 + (X22*len)/N + ((X31)/(N^4) + (X32)/(N^2) + X33)*(4*len^2*(tm - 2*t1/3 - 2*t3/3)/(2*tm - t1 - t3)^2) + ((X41)/(N^4) + (X42)/(N^2) + X43)*(4*len^2*(t1+t3)/((2*tm - t1 - t3)^2 * t1 * t3));


N = (1:0.1:50);
Es = (X11*tm)./N.^2 + (X21*len)./N.^3 + (X22*len)./N + ((X31)./(N.^4) + (X32)./(N.^2) + X33)*(4*len^2*(tm - 2*t1/3 - 2*t3/3)/(2*tm - t1 - t3)^2) + ((X41)./(N.^4) + (X42)./(N.^2) + X43)*(4*len^2*(t1+t3)/((2*tm - t1 - t3)^2 * t1 * t3));
plot(N, Es, 'b-');grid on; title('E against N');