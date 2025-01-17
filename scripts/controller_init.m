X1 = R*D^2/Ke^2;
X2 = R*2*B*D/Ke^2 + D;
X3 = R*B^2/Ke^2 + B;
X4 = (R*A^2 + A*B*L)/Ke^2;

% trapezoidal, cosine, cubic
MP = 'trapezoidal';

% len = 30;
% tm = 8;
% k = 0.3552;
% 
% max_vel = 8.9;
% max_acc = 10;

len = 30;
tm = 10;
k = 0;

max_vel = 8.9;
max_acc = 10;


% --------------------------- ONLY ONE SETTING SHOULD BE ON
USE_ACC_K = 0;
USE_ACC_TM = 0;
USE_ACC = 0;

USE_VEL_K = 0;
USE_VEL_TM = 0;
USE_VEL = 0;

USE_ACC_VEL = 0;

USE_K = 0;
USE_TM = 0;

USE_K_TM = 1;
% --------------------------- ONLY ONE SETTING SHOULD BE ON


mp = 0;
if(strcmp(MP,'trapezoidal') == 1)
    mp = 1;
    if(USE_ACC_K+USE_ACC_TM+USE_ACC+USE_VEL_K+USE_VEL_TM+USE_VEL+USE_ACC_VEL+USE_K+USE_TM+USE_K_TM ~= 1)
        error('INVALID SETTING');
    end

    if(USE_ACC_K)
        tm = sqrt(4*len/(max_acc*(1-k^2))); %#ok<*UNRCH>
        max_vel = (2/1+k) * (len/tm);
    elseif(USE_ACC_TM)
        k = sqrt(1-(4*len/(max_acc*tm^2)));
        max_vel = (2/1+k) * (len/tm);
        % elseif(USE_ACC)
    elseif(USE_VEL_K)
        tm = 2*len/(max_vel*(1+k));
        max_acc = max_vel/(tm*(1-k)/2);
    elseif(USE_VEL_TM)
        k = (2*len/(max_vel*tm))-1;  % (FIXED)Wrong should not take abs, the value should be positive below some absolute max vel....trapz eqns allow for only same accel and decel at opt
        if(k<0)
            error('MAX_VEL too high to meet the motion profile requirments');
        elseif(k>=1)
            error('MAX_VEL too low to meet the motion profile requirements');
        end
        max_acc = max_vel/(tm*(1-k)/2);
        % elseif(USE_VEL)
    elseif(USE_ACC_VEL)
        k = roots([(max_vel^2+max_acc*len) (2*max_vel^2) (max_vel^2-max_acc*len)]);
        k = k(k<1&k>=0);
        if(isempty(k))
            error('Unsuitable MAX_VEL and MAX_ACC values');
        end
        tm = 2*len/(max_vel*(1+k));
    elseif(USE_K)
        tm = roots([X1*(3*(1-k)/2 - 9*(1-k)^2/4 + 9*(1-k)^3/8 - 3*(1-k)^4/16) (0) X3*len^2*(5*(1-k)^2/4 - 3*(1-k)/2) (0) (12*X4*len^2)]);
        tm = abs(tm(1));
        max_vel = (2/1+k) * (len/tm);
        max_acc = (4/(1-k^2)) * (len/tm^2);
    elseif(USE_TM)
        k = roots([X3*tm^2 -2*X3*tm^2 ((X3*tm^2)-(18*X4)) 6*X4]);
        k = min(k(k<1&k>=0));
        max_vel = (2/1+k) * (len/tm);
        max_acc = (4/(1-k^2)) * (len/tm^2);
    elseif(USE_K_TM)
        % normal case
        max_vel = (2/(1+k)) * (len/tm);
        max_acc = (4/(1-k^2)) * (len/tm^2);
    end
elseif(strcmp(MP,'cosine') == 1)
    mp = 2;
    if(USE_ACC_VEL)
        cof = (2*pi)/(max_acc*4*len);
        k = roots([(cof*max_vel^2+1) (2*cof*max_vel^2) (cof*max_vel^2-1)]);
        k = k(k<1&k>=0);
        if(isempty(k))
            error('Unsuitable MAX_VEL and MAX_ACC values');
        end
        tm = 2*len/(max_vel*(1+k));
    elseif(USE_K_TM)
        % normal case
        max_acc = (2*pi/(1-k^2))*(len/tm^2);
        max_vel = (2/(1+k)) * (len/tm);
    end
elseif(strcmp(MP,'cubic') == 1)
    mp = 3;
    if(USE_ACC_VEL)
        cof = (-12*max_vel^2)/(max_acc*9*len);
        k = roots([(cof-1) (4*cof-1) (4*cof+2)]);
        k = k(k<1&k>=0);
        if(isempty(k))
            error('Unsuitable MAX_VEL and MAX_ACC values');
        end
        tm = 3*len/(max_vel*(2+k));
    elseif(USE_K_TM)
        % normal case
        max_acc = (-12/(k^2+k-2))*(len/tm^2);
        max_vel = (3/(2+k)) * (len/tm);
    end
end


ctrl_max_acc = max_acc;
ctrl_max_vel = max_vel;

ctrl_len = len;
ctrl_tm = tm;
ctrl_k = k;

ctrl_t1 = (tm-k*tm)/2;
ctrl_t3 = ctrl_t1;

ctrl_mp = mp;



