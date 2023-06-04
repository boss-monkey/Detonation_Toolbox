%% 斜爆轰发动机理论推力计算
%author :Boss Monkey
%Email : Baiht0201@nuaa.edu.cn 
%参考文献 斜爆轰推进理论、技术及其实验验证
clear all
%% 参数输入
load('data_nasa9.mat')
load('.\molecular_w.mat') 
m = ["H2", "O2" , "N2", "H2O" , "O" , "H", "OH", "H2O2","HO2","N","NO"];
c_air = [0, 0.21, 0.79, 0, 0, 0, 0, 0, 0, 0, 0]; %来流摩尔数
c_flue = [0.42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];   %燃料摩尔数
A = [2, 0, 0, 2, 0, 1, 1, 2, 1, 0, 0; %氢原子数
     0, 2, 0, 1, 1, 0, 1, 2, 2, 0, 1; %氧原子数
     0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 1]; %氮原子数
Ru = 8.314;
max = 100;
e =10^(-7);
tic
%% 计算基本量
ns = length(m);
Mw = zeros(1,ns);   %分子摩尔量
for i = 1:ns
    Mw(i)=m_w.(m(i));
end
D = zeros(ns,9,3);     %Extract NASA 9 coefficients
for j = 1:ns
    D(j,:,:) = (coeff_nasa9.(m(j)))';
end
h = zeros(1,21);
v = h;
for k = 1:21
%% 斜激波求解
Ma = 8+(k-1)*7/20;                     %飞行马赫数（Ma>7）               
p_0 = 1197;                  %大气静压
T_0 = 226.51;                %大气静温
theta_b = deg2rad(7.5);       %进气道压缩角
Mm_0 = dot(c_air,Mw)/sum(c_air); 
Rm_0 = Ru/Mm_0;
rho_0 = p_0/(Rm_0*T_0);
a_s_0 = a_s (T_0,Mm_0,Rm_0,c_air,D);
u_0 = Ma * a_s_0; 
A_0 = Mm_0*sum(c_air) /(rho_0*u_0);   %进气捕获面积
%% 波后参数
Mm_1 = Mm_0;
Rm_1 = Rm_0;

if theta_b == 0
    T_1 = T_0;
    p_1 = p_0;
    beta_0 =0;
    rho_1 = rho_0;
    u_1 = u_0;
else 
    outlet = shock_o(u_0, p_0 ,T_0, theta_b, Mm_0, Mm_0, c_air ,c_air, D);
    T_1 = outlet(1);
    p_1 = outlet(2);
    beta_0 = rad2deg(outlet(3));
    rho_1 = p_1/(Rm_1*T_1);
    u_1 = (u_0*cos(outlet(3))*tan(outlet(3)-theta_b))/sin(outlet(3)-theta_b);
end 

A_1 = Mm_1*sum(c_air)/(rho_1*u_1); %进气道出口面积
a_s_1 = a_s (T_1,Mm_1,Rm_1,c_air, D);
Ma_1 = u_1/a_s_1;
%% 燃料混合
flue_T_t = 600;          %燃料总温
flue_Ma = 2;             %喷射马赫数
Af = 0.01;               %燃料面积与流道的面积比,人为指定

flue_cal = total_cal (c_flue, flue_Ma, flue_T_t, Mw, D);
phi = 1;                 %当量比
c = c_air+phi*c_flue;    %混合后气体成分

flue_Mm = flue_cal(1); 
flue_Rm = flue_cal(2);
flue_T = flue_cal(5);

flue_a_s = a_s (flue_T,flue_Mm,flue_Rm,c_flue,D);
flue_u = flue_Ma * flue_a_s;
mf = phi*flue_Mm*sum(c_flue)/(Mm_1*sum(c_air));   %质量比
flue_rho = rho_1 * u_1 * mf / (flue_u*Af);
flue_p = flue_rho*flue_Rm*flue_T;                 %喷射出口压强
flue_p_total = flue_p/flue_cal(4);                %燃料总压

Mm_2 = dot(c,Mw)/sum(c);                %混合后气体相对质量
Rm_2 = Ru / Mm_2;

Z_M_0 = (rho_1*u_1+flue_rho*flue_u*Af)/(Af+1);
Z_F_0 = ((rho_1*u_1^2 +p_1)+(flue_rho*flue_u^2+flue_p)*Af)/(Af+1);
Z_H_0 = ((h0_m(T_1,c_air,Mm_1,D) +u_1^2 /2)*rho_1*u_1+...
      (h0_m(flue_T,c_flue,flue_Mm,D) +flue_u^2 /2)*...
      flue_rho*flue_u*Af)/(rho_1*u_1+flue_rho*flue_u*Af);

u = u_1;
for i=1:100
    T_2 = 1/Rm_2 * (Z_F_0/Z_M_0 - u)*u;
    if T_2 <200
        T_2 = 200;
    end
    f = h0_m(T_2,c,Mm_2,D) +u^2 /2 - Z_H_0;
    g = Cp_m(T_2,c,Mm_2,D)/Rm_2*(Z_F_0/Z_M_0 - 2*u)+u;
    u_2 = u - f/g;
    if abs(u-u_2)<e
        break
    end
    u=u_2;
end
T_2;
rho_2 = (rho_1*u_1+flue_rho*flue_u*Af)/((Af+1)*u_2);
p_2 = rho_2*Rm_2*T_2;
A_2 = Mm_2*sum(c)/(rho_2*u_2);      %混合后流道面积
%% 反应激波求解
theta_w = deg2rad(18);              %爆轰楔面角

s_0 = shock_o(u_2, p_2 ,T_2, theta_w, Mm_2, Mm_2, c, c, D);
f_0 = post_reaction(s_0(1),s_0(1), c, D, A);

 for i = 1:max
     Mm_3 = dot(f_0,Mw)/sum(f_0);
     s_1 = shock_o(u_2, p_2 ,T_2, theta_w, Mm_2, Mm_3, c,f_0, D);
     f_1 = post_reaction(s_1(1), s_1(2),c, D, A);
    if abs(f_1(1)-f_0(1))< e
        break
    end
    f_0 = f_1 + abs(f_1-f_0)*0.7;
 end

%% 波后参数 
beta_1 = rad2deg(s_1(3));           %爆轰波波角
Rm_3 = Ru/Mm_3;
T_3 = s_1(1);
p_3 = s_1(2);
rho_3 = p_3/(Rm_3 * T_3); 
u_3 = (u_2*cos(s_1(3))*tan(s_1(3)-theta_w))/sin(s_1(3)-theta_w);
c_3 = f_0;                  %波后物质
A_3 = Mm_3*sum(c_3)/(rho_3*u_3);   %燃烧室出口面积
%% 等熵膨胀求解
A_t = A_0/A_3;                      %喷管扩张面积比，人为给定
e_0 = isentropic_cal_area(p_3, A_t, T_3, u_3, c_3, c_3, Mw, D);  %指定喷管面积
% e_0 = isentropic_cal_pressure(p_3, p_0, T_3, u_3, c_3, c_3, Mw, D);  %指定出口压强
g_0 = post_reaction(e_0(1),e_0(2), c, D, A );
 for i = 1:max
    e_1 = isentropic_cal_area(p_3, A_t, T_3, u_3, c_3, g_0, Mw, D);
    %e_1 = isentropic_cal_pressure(p_3, p_0, T_3, u_3, c_3, g_0, Mw, D);
    g_1 = post_reaction(e_1(1), e_1(2), c, D, A);
    if abs(g_1(1)-g_0(1))< e
        break
    end
    g_0 = g_1 + abs(g_1-g_0)*0.4;
 end
%% 喷管后流场参数 
T_4 = e_1(1);
p_4 = e_1(2);
u_4 = e_1(3);
c_4 = g_0;
Mm_4 = dot(g_0,Mw)/sum(g_0); 
rho_4 = p_4/(Ru/Mm_4 * T_4);
A_4 = Mm_4*sum(c_4)/(rho_4*u_4);   %喷管出口面积（此时应等于进气捕获面积）
%% 发动机性能
T_e = (rho_4*(cos(theta_w-theta_b)*u_4)^2+p_4)*A_4 - (rho_0*u_0^2+p_0)*A_0;  %发动机单位推力
m_a = rho_0*u_0*A_0;                                  %空气质量流量
m_f = flue_rho*flue_u*A_2/(Af+1)*Af;                  %燃料质量流量
I_s = T_e/(m_f*9.8);                                  %比冲
h(k)=I_s;
q(k)=T_4;
v(k)=Ma;
end

plot(v,h);

toc

