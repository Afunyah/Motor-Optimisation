if(~NONE)
    [OPTXS,FVAL,~,OUTPUT] = fmincon(funE,OPTXS,[],[],[],[],lb,ub,nonlinfcn,hybridopts);
    cmes2 = OUTPUT.message;

    tm = OPTXS(1); t1 = OPTXS(2); t3 = OPTXS(3); max_vel = OPTXS(4);
    k = (tm-t1-t3)/tm;
    E_total = FVAL;
end


calc_opt_res = sprintf('tm = %0.2f s \nt1 = %0.2f s \nt3 = %0.2f s \nmax_vel = %0.2f rad/s \nk = %0.2f\nE = %0.3f J',tm,t1,t3,max_vel,k,E_total);
calc_opt_message = cat(2,cmes1,cmes2,calc_opt_res);

run("controller_calc.m");