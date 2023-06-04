%燃料喷射参数的转换,滞止参数与静参数的转换
function total_cal = total_cal (c, Ma, T_t, Mw, D )
Ru = 8.314;
Mm = dot(c,Mw)/sum(c);
Rm = Ru/Mm;
gamma = Cp_m(T_t,c,Mm,D)/(Cp_m(T_t,c,Mm,D)-Rm);
%P_0/P
pi = (1 + (gamma-1) /2 * Ma^2 ) ^((-1) * gamma / (gamma-1));
%T_0/T
tau = (1 + (gamma-1) /2 * Ma^2 )^(-1);   
T = T_t* tau;
total_cal =[Mm,Rm,tau,pi,T];
end
