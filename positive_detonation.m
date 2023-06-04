%过驱正爆轰波的计算
%author :Boss Monkey
%Email : Baiht0201@nuaa.edu.cn 
clear all
load('data_nasa9.mat')
load('.\molecular_w.mat') ;
%% 参数输入
c = [0.42, 0.21, 0.79, 0, 0 ,0, 0, 0, 0, 0, 0 ];
m = ["H2", "O2" , "N2", "H2O" , "O" , "H", "OH", "H2O2","HO2","N","NO"];
A = [2, 0, 0, 2, 0, 1, 1, 2, 1, 0, 0; %氢原子数
     0, 2, 0, 1, 1, 0, 1, 2, 2, 0, 1; %氧原子数
     0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 1]; %氮原子数
T_1 = 300; p_1 = 101325; u_1 = 2500;
%% 基本量计算
ns = length(m);
Ru = 8.314;
Mw = zeros(1,ns);  
for i = 1:ns
    Mw(i)=m_w.(m(i));
end
D = zeros(ns,9,3);     %Extract NASA 9 coefficients
for j = 1:ns
    D(j,:,:) = (coeff_nasa9.(m(j)))';
end
Mm_1 = dot(c,Mw)/sum(c); 
Rm_1 = Ru/Mm_1;
a_0 = shock_n(u_1, p_1 ,T_1, Mm_1, Mm_1, c, c, D);
b_0 = post_reaction(a_0(1),a_0(2), c, D, A);
max = 100;
e =10^(-6);
%% 激波求解
 for i = 1:max
    Mm_2 = dot(b_0,Mw)/sum(b_0); 
    a_1 = shock_n(u_1, p_1 ,T_1, Mm_1, Mm_2,c, b_0, D);
    b_1 = post_reaction(a_1(1), a_1(2), c, D, A);
    if abs(b_1(1)-b_0(1))< e
        break
    end
    b_0 = b_1 + abs(b_1-b_0)*0.4;
 end
%% 输出 
rho_1 = p_1/(Rm_1*T_1);
T_2 = a_1(1);
p_2 = a_1(2);
u_2 = a_1(3);
Rm_2 = Ru/Mm_2;
c_2 = b_1;
rho_2 = p_2/(Rm_2 * T_2); 

