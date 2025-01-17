
syms S P M H G Q C W POS1 POS2 POS3 VEL1 VEL2 VEL3 ACC1 ACC2 ACC3 JERK1 JERK2 JERK3
syms t ta tb tc tm t1 t3 A1 A2
syms A B D R L Ke len

% JERK1 = 0;
% JERK2 = 0;
% JERK3 = 0;
% 
% ACC1 = A1;
% ACC2 = 0;
% ACC3 = -A2;
% 
% VEL1 = A1 * t;
% VEL2 = A1 * tb;
% VEL3 = (A1 * tb) - A2*(t-tc);
% 
% POS1 = 0.5 * A1 * t^2;
% POS2 = (0.5 * A1 * tb^2) + (A1 * tb) * (t - tb);
% POS3 = (0.5 * A1 * tb^2) + (A1 * tb) * (tc - tb) + (A1 * tb) * (t - tc) + (A2 * tc) * (t - tc) - 0.5 * A2 * (t^2 - tc^2);

POS1 = 0.5 * A1 * t^2;
POS2 = (0.5 * A1 * tb^2) + (A1 * tb) * (t - tb);
POS3 = (0.5 * A1 * tb^2) + (A1 * tb) * (tc - tb) + (A1 * tb) * (t - tc) + (A2 * tc) * (t - tc) - 0.5 * A2 * (t^2 - tc^2);

VEL1 = diff(POS1, t, 1);
VEL2 = diff(POS2, t, 1);
VEL3 = diff(POS3, t, 1);

ACC1 = diff(POS1, t, 2);
ACC2 = diff(POS2, t, 2);
ACC3 = diff(POS3, t, 2);

JERK1 = diff(POS1, t, 3);
JERK2 = diff(POS2, t, 3);
JERK3 = diff(POS3, t, 3);



P1 = S + P*VEL1 + M*ACC1 + H*VEL1*ACC1 + G*VEL1^2 + Q*ACC1^2 + C*VEL1*JERK1 + W*JERK1;
P2 = S + P*VEL2 + M*ACC2 + H*VEL2*ACC2 + G*VEL2^2 + Q*ACC2^2 + C*VEL2*JERK2 + W*JERK2;
P3 = S + P*VEL3 + M*ACC3 + H*VEL3*ACC3 + G*VEL3^2 + Q*ACC3^2 + C*VEL3*JERK3 + W*JERK3;

E1 = int(P1,t,[ta tb]);
E2 = int(P2,t,[tb tc]);
E3 = int(P3,t,[tc tm]);

E4 = E1+E2+E3;

E5 = subs(E4, {tb, tc}, {t1-ta,tm-t3});

S = R*D^2/Ke^2;
P = 2*B*D*R/Ke^2 + D;
M = (2*A*D*R + D*B*L)/Ke^2;
H = (2*A*B*R + B^2*L)/Ke^2 + A;
G = B^2*R/Ke^2 + B;
Q = (A^2*R + A*B*L)/Ke^2;
C = (A^2+A*B)*L/Ke^2;
W = (A*D)*L/Ke^2;

E6 = subs(E5);
E7 = collect(E6);
E8 = simplify(E7);
E9 = collect(E8,[S,P,M,H,G,Q,C,W]);

E = E9;

% len = (2*tm-t1-t3)*A1*t1/2;

