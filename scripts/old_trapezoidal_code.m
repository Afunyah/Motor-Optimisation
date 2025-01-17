% 
% len = 5*pi;
% tm = 2;
% 
% 
% % max_ang_acc = 0.01;
% % k = sqrt(1 - (4*len/tm^2 * max_ang_acc));
% 
% k = 0.5;
% 
% max_ang_vel = (2/(1+k)) * (len/tm);
% % max_ang_vel = 0; k = (2*L/(tm*max_ang_vel)) - 1;
% 
% ta = 0;
% tb = ta + 0.5 * (tm - k*tm);
% tc = tb + k*tm;
% td = tm;
% % td = tc + tb;
% 
% alpha = (4/(1-k^2)) * (len/tm^2);
% 
% T = ta:0.01:tm;
% ang_acc_arr = zeros(size(T));
% ang_vel_arr = zeros(size(T));
% ang_pos_arr = zeros(size(T));
% 
% for i=1:length(T)
%     t = T(i);
%     if (t>=ta) && (t<=tb)
%         ang_acc = alpha;
%         ang_vel = alpha * t;
%         ang_pos = 0.5 * alpha * t^2; 
%     elseif (t>tb) && (t<=tc)
%         ang_acc = 0;
%         ang_vel = alpha * tb;
%         ang_pos = (0.5 * alpha * tb^2) + (alpha * tb) * (t - tb); 
%     elseif (t>tc) && (t<=td)
%         ang_acc = -alpha;
%         ang_vel = (alpha * tb) - alpha*(t-tc);
%         ang_pos = (0.5 * alpha * tb^2) + (alpha * tb) * (tc - tb) + (alpha * tb) * (t - tc) + (alpha * tc) * (t - tc) - 0.5 * alpha * (t^2 - tc^2);
%     end
%     
%     ang_acc_arr(i) = ang_acc;
%     ang_vel_arr(i) = ang_vel;
%     ang_pos_arr(i) = ang_pos;
% 
%     
% end
% 
% 
% plot(T,ang_vel_arr,'b-');
% hold on;
% plot(T,ang_acc_arr,'r-');
% plot(T,ang_pos_arr,'g-');
% 
% 
% 
% 
% 

function [acc,vel,pos] = fcn(t)

    len = 5*pi;
    tm = 4;
    
    
    % max_ang_acc = 0.01;
    % k = sqrt(1 - (4*len/tm^2 * max_ang_acc));
    
    k = 0.5;
    
%     max_ang_vel = (2/(1+k)) * (len/tm);
    % max_ang_vel = 0; k = (2*L/(tm*max_ang_vel)) - 1;
    
    ta = 0;
    tb = ta + 0.5 * (tm - k*tm);
    tc = tb + k*tm;
    td = tm;
    
    alpha = (4/(1-k^2)) * (len/tm^2);
    
    if (t>=ta) && (t<=tb)
        ang_acc = alpha;
        ang_vel = alpha * t;
        ang_pos = 0.5 * alpha * t^2; 
    elseif (t>tb) && (t<=tc)
        ang_acc = 0;
        ang_vel = alpha * tb;
        ang_pos = (0.5 * alpha * tb^2) + (alpha * tb) * (t - tb); 
    elseif (t>tc) && (t<=td)
        ang_acc = -alpha;
        ang_vel = (alpha * tb) - alpha*(t-tc);
        ang_pos = (0.5 * alpha * tb^2) + (alpha * tb) * (tc - tb) + (alpha * tb) * (t - tc) + (alpha * tc) * (t - tc) - 0.5 * alpha * (t^2 - tc^2);
    else
        ang_acc = 0;
        ang_vel = 0;
        ang_pos = (0.5 * alpha * tb^2) + (alpha * tb) * (tc - tb) + (alpha * tb) * (tm - tc) + (alpha * tc) * (tm - tc) - 0.5 * alpha * (tm^2 - tc^2);
    end
    


acc = ang_acc;
vel = ang_vel;
pos = ang_pos;

