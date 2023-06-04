%alpha-QSS方法计算单点上的化学反应
function cal = cal_Qc(t, rho, p, y_i, reaction,mix) 
Ru = 8.31442;
Rc = 1.98718;
Ra = 82.057338;  %和大气压相关的普适气体常数；

v_1 = reaction.v_1; v_2 = reaction.v_2; A =reaction.A; beta = reaction.beta;
Ea = reaction.Ea; Mr = reaction.Mr; nr =reaction.nr; nt =reaction.nt;
Mw =mix.Mw; ns= mix.ns;D = mix.D;

e = ones(1,ns)*10^-16;
e1 = ones(1,ns)*10^-3;
scale = 10^-6;
one = ones(1,ns);
%% Calculate the basic quantity
R_i = Ru./ Mw;
Rm =  sum(R_i.*y_i);
T = p/(Rm * rho);

data_cal = data_cal_vector(T,D);
cp_i = data_cal(1,:)./Mw;
h_i = data_cal(2,:)./Mw;
s_i = data_cal(3,:)./Mw;

T0 = T;
y_i0 = y_i;
%% Reaction
kf = A .* T.^beta .* exp(-Ea/(Rc*T));

c_i = (rho*y_i./Mw* scale);
F = c_i.^v_1;
B = c_i.^v_2;

Gr = ones(nr,1);
Gr(1:nt) =sum(Mr.*c_i,2);

vr= v_2 - v_1;
kpr = exp(sum((vr.*(s_i./R_i)),2) - sum((vr.* h_i./(R_i*T)),2));
kcr = (1/Ra/T).^sum(vr,2).*kpr;     
kb = kf./kcr;     %Equilibrium constant

%% THE alpha-QSS ALGORITHM
flag = 1;
dt = 1e-8;
t_final = 0;
if t<dt
    dt = t;
end
while(flag)
    qi =  1/scale * Mw .* sum(Gr.* (v_2.* kf .* prod(F,2)+ v_1 .* kb .* prod(B,2)))/rho;
    pi =  1./(c_i+e) .* sum(Gr.* (v_1.* kf .* prod(F,2)+ v_2 .* kb .* prod(B,2)));

    alpha_i = (180* one + 60*(pi*dt) + 11*(pi*dt).^2 + (pi*dt).^3)./...
            (360* one + 60*(pi*dt) + 12*(pi*dt).^2 + (pi*dt).^3);

    y0 = y_i; q0 = qi; p0 = pi; a0 = alpha_i;
    % Predictor
    yp = y0 + (q0-p0.*y0)*dt/(one +a0.*p0*dt);

    %Calculate the parameter for the corrector step
    c_ip = (rho*yp./Mw* scale);
    F = c_ip.^v_1;
    B = c_ip.^v_2;

    Gr(1:nt) =sum(Mr.*c_ip,2);
    q_p =  1/scale * Mw .* sum(Gr.* (v_2.* kf .* prod(F,2)+ v_1 .* kb .* prod(B,2)))/rho;
    p_p =  1./(c_ip+e) .* sum(Gr.* (v_1.* kf .* prod(F,2)+ v_2 .* kb .* prod(B,2)));

    p_c = (p0+p_p)/2;
    a_p = (180* one + 60*(p_c*dt) + 11*(p_c*dt).^2 + (p_c*dt).^3)./...
          (360* one + 60*(p_c*dt) + 12*(p_c*dt).^2 + (p_c*dt).^3);
    q_c = q_p .* a_p + (one - a_p).* q0;
    %Corrector
    yc = y0 + (q_c-p_c.*y0)*dt./(one +a_p.*p_c*dt);
    

    %Timestep Selection
    c0 = 1.6;   %constant to control precision
    sigma = max(abs((yc - yp)./((yc+e1)*c0*10^-3)));     %Attention! The parameters here will affect the stabilization

    if sigma < 2
        dt = dt*min([sigma^-0.5,1.6]);      % constant = 1.6 to limit the step size
        if t_final+ dt >t
            dt = t - t_final;
            flag = 0;
        end
        t_final = t_final+ dt ;
        %Calculate the temperature for next step
        dydt  = (q_c-p_c .* y0)./(one +a_p.*p_c*dt);
        dTdt = (sum(h_i.*dydt)-T*sum(R_i.*dydt))/(-sum(cp_i.*y0)+sum(R_i.*y0));
        T = T+ dTdt*dt;

        %Calculate the reactions again
        kf = A .* T.^beta .* exp(-Ea/(Rc*T));
        y_i = yc;
        c_i = (rho*yc./Mw* scale);
        F = c_i.^v_1;
        B = c_i.^v_2;
        Gr(1:nt) =sum(Mr.*c_i,2);

        data_cal = data_cal_vector(T,D);
        cp_i = data_cal(1,:)./Mw;
        h_i = data_cal(2,:)./Mw;
        s_i = data_cal(3,:)./Mw;
       
        %Equilibrium constant
        kpr = exp(sum((vr.*(s_i./R_i)),2) - sum((vr.* h_i./(R_i*T)),2));
        kcr = (1/Ra/T).^sum(vr,2).*kpr;     
        kb = kf./kcr;
        

        
    else
        dt = dt/3;
    end
    if  dt < 1e-16   %Avoid getting stuck
        T = T0;
        y_i = y_i0;  
        flag = 0;
    end
  
end
    
    gamma = cal_gamma_T(T,y_i,mix);
    R_1 =  sum(R_i.*y_i);
    
cal = [rho*R_1*T,gamma,T,y_i];   %（p，T，gamma，y）
