function shock_n = shock_n(u_1, p_1, T_1 ,Mm_1, Mm_2, c, y, D )
max = 20;

x_0 = 400;%起始点
e=10^(-6);%误差范围
R_u =8.314;
Rm_1 = R_u/Mm_1;
Rm_2 = R_u/Mm_2;
%% 计算
for i=1:max
    T = 1 / Rm_2 *(Rm_1 *T_1 /u_1 +u_1 -x_0)*x_0;
    f = h0_m(T,y,Mm_2,D) +x_0^2 /2 - h0_m(T_1,c,Mm_1,D) - u_1^2 /2;
    g = Cp_m(T,y,Mm_2,D)/Rm_2 *(Rm_1*T_1/u_1 + u_1 - 2*x_0)+x_0;
    x_1= x_0-f/g;
    if abs(x_1-x_0)<e
        break
    end
    x_0=x_1;
end

%% 输出
rho_1 = p_1/(Rm_1*T_1);
T_2 = T;
u_2 = x_0;
rho_2 = rho_1 *u_1/u_2;
p_2 = rho_2 * Rm_2 * T_2;
shock_n = [T_2, p_2, u_2];