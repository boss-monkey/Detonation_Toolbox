function dydt =ode_reaction(t, y, reaction,mix, rho) 

Ru = 8.31442;
Rc = 1.98718;
Ra = 82.057338;  %和大气压相关的普适气体常数；
pa = 101325;
v_1 = reaction.v_1; v_2 = reaction.v_2; A =reaction.A; beta = reaction.beta;
Ea = reaction.Ea; Mr = reaction.Mr; nr =reaction.nr; nt =reaction.nt;
Mw =mix.Mw; ns= mix.ns;D = mix.D;
R_i = Ru./ Mw;


data_cal = data_cal_vector(y(1),D);
cp_i = data_cal(1,:)./Mw;
h_i = data_cal(2,:)./Mw;
s_i = data_cal(3,:)./Mw;


R = dot(R_i,y(2:end)');
p = rho*R*y(1); 

kf = A .*  y(1).^beta.* exp(-Ea/(Rc*y(1)));

c_i = (rho* y(2:end)'./Mw * 10^-6);
F = c_i.^v_1;
B = c_i.^v_2;

Gr = ones(nr,1);
Gr(1:nt) =sum(Mr.*c_i,2);

vr= (v_2 - v_1);
kpr = exp(sum((vr.*(s_i./R_i)),2) - sum((vr.* h_i./(R_i*y(1))),2));
%kcr = (1/Ra/y(1)).^sum(vr,2) .*kpr;  
kcr = (p/pa/sum(c_i)).^(-sum(vr,2)).*kpr;
kb = kf./kcr;

   
o = 10^6* Mw.* Gr.* vr.* (kf .* prod(F,2)- kb .* prod(B,2));  % 各反应各组分生成速率
o_i = sum(o);

dy_idt = 1/rho*o_i;
dTdt = (dot(h_i,dy_idt)-y(1)*dot(R_i,dy_idt))/(-dot(cp_i,y(2:end)')+dot(R_i,y(2:end)'));
dydt = [dTdt;zeros(ns,1)];
for i = 1:ns
    dydt(i+1) = dy_idt(i);
end

