%
function dydx = ode_ZND(x, y, reaction,mix,flux, D) 

Ru = 8.31442;
Rc = 1.98718;
Ra = 82.057338;  %和大气压相关的普适气体常数；
pa = 101325;
v_1 = reaction.v_1; v_2 = reaction.v_2; A =reaction.A; beta = reaction.beta;
Ea = reaction.Ea; Mr = reaction.Mr; nr =reaction.nr; nt =reaction.nt;
 Mw =mix.Mw; ns= mix.ns;
ma = flux(1);mo =flux(2); 

rho = ma/y(1);
p = mo - rho*y(1)^2;
R_i = Ru./ Mw;
R = dot(R_i,y(2:end)');
T = p /rho/R;


data_cal = data_cal_vector(T,D);
cp_i = data_cal(1,:)./Mw;
h_i = data_cal(2,:)./Mw;
s_i = data_cal(3,:)./Mw;
cp = dot(cp_i,y(2:end)');


kf = A .* T.^beta .* exp(-Ea/(Rc*T));

c_i = (rho* y(2:end)'./Mw * 10^-6);
F = c_i.^v_1;
B = c_i.^v_2;

Gr = ones(nr,1);
Gr(1:nt) =sum(Mr.*c_i,2);

vr= (v_2 - v_1);
kpr = exp(sum((vr.*(s_i./R_i - h_i./(R_i*T))),2));
kcr = (1/Ra/T).^sum(vr,2).*kpr;  
%kcr = (p/pa/sum(c_i)).^(-sum(vr,2)).*kpr;
kb = kf./kcr;
sum(c_i)   
o = 10^6* Mw.* Gr.* vr.* (kf .* prod(F,2)- kb .* prod(B,2));  % 各反应各组分生成速率
o_i = sum(o);

dy_idx = 1/ma *o_i;
dudx = (dot(h_i,dy_idx) - cp*T/R * dot(R_i,dy_idx))/(cp*y(1)/R - cp*T/y(1)-y(1));

dydx = [dudx;zeros(ns,1)];
for i = 1:ns
    dydx(i+1) = dy_idx(i);
end

