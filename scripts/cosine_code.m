function [acc,vel,pos] = cosine_code(t)

    len = 5*pi;
    tm = 4;
        
    k = 0.5;
    
    ta = 0;
    tb = ta + 0.5 * (tm - k*tm);
    tc = tb + k*tm;
    td = tm;

    t_acc = tb - ta;
    t_dec = td - tc;
    
%     max_ang_acc = (2*pi/(1-k^2))*(len/tm^2);
    max_ang_vel = (2/(1+k)) * (len/tm);

    thet1 = pi/t_acc;
    thet3 = pi/t_dec;
    
    if (t>=ta) && (t<=tb)
        ang_acc = max_ang_vel*(thet1*sin(thet1*t))/2;
        ang_vel = max_ang_vel*(1-cos(thet1*t))/2;
        ang_pos = max_ang_vel*(t-sin(thet1*t)/thet1)/2; 
    elseif (t>tb) && (t<=tc)
        ang_acc = 0;
        ang_vel = max_ang_vel;
        ang_pos = max_ang_vel*(t-t_acc/2); 
    elseif (t>tc) && (t<=td)
        tprm = t - tm + t_dec;
        ang_acc = -max_ang_vel*(thet3*sin(thet3*tprm))/2;
        ang_vel = max_ang_vel*(1+cos(thet3*tprm))/2;
        ang_pos = max_ang_vel*(tm - t_acc/2 - t_dec) + max_ang_vel*(tprm+sin(thet3*tprm)/thet3)/2;
    else
        tprm = tm - tm + t_dec;
        ang_acc = 0;
        ang_vel = 0;
        ang_pos = max_ang_vel*(tm - t_acc/2 - t_dec) + max_ang_vel*(tprm+sin(thet3*tprm)/thet3)/2;
    end
    


acc = ang_acc;
vel = ang_vel;
pos = ang_pos;

