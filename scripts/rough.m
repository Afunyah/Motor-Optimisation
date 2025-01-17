% % ppp = out(1,1:end).motor_data.E.Data(end);
% ccc = size(out);
% for i=1:length(out)
%     ccc(i) = out(1,i).motor_data.E.Data(end);
% end
% 
% % [xx,yy,zz] = meshgrid((1:0.3:10),(0.1:0.02:0.8),ccc);
% 
% zz = reshape(ccc, [36 31]);
% [xx,yy] = meshgrid((1:0.3:10),(0.1:0.02:0.8));

% syms x
% expr = x*log(1+x);
% F = int(expr,[0 1])

syms t ta tb T A1 A2 S B D F G

F0 = int(S,t,[0 T]);

F1 = int(B*A1*t + D*A1 + F*A1^2*t^2 + G*A1^2,t,[0 ta]);

F2 = int(B*A1*ta + F*A1^2*ta^2,t,[ta tb]);

F3 = int(B*A1*ta - B*A2*t + B*A2*tb - D*A2 + F*A1^2*ta^2 + F*A2^2*t^2 + F*A2^2*tb^2 - 2*F*A1*A2*ta*t + 2*F*A1*A2*ta*tb - 2*F*A2^2*tb*t + G*A2^2,t,[tb T]);

Fint = F0+F1+F2+F3;

syms t1 t3 
Fint2 = subs(Fint, T-tb, t3);
Fint2 = subs(Fint2, tb, T-t3);
Fint2 = subs(Fint2, ta, t1);
Fint3 = collect(Fint2, [B,D,F,G]);
% The form above is good, for using multiple variables
% Remember to add gear ratios and so on if i want to add them

% Doing verification tests with the eqns below
syms L
Fint4 = subs(Fint2, A2, A1);
Fint5 = subs(Fint4, A1, 2*L/((2*T-t1-t3)*t3));
Fint6 = collect(Fint5, [B,D,F,G]);
Fint7 = collect(Fint6);
Fint8 = collect(Fint7, [S,B,D,F,G]);

syms X1 X2 X3 X4 tm len
Fint9 = subs(Fint8,{S,B,F,G,D},{X1,X2,X3,X4,0});
Fint10 = collect(Fint9,[X1,X2,X3,X4]);
Fint11 = collect(Fint10);

% t1 = t3
% Fint12 = subs(Fint11, t3,t1);

Fint13 = subs(Fint11, {T,L},{tm,len});



