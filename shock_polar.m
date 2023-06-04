%斜爆轰波驻定窗口的计算
%author :Boss Monkey
%Email : Baiht0201@nuaa.edu.cn 
%脱体点附近会发散，采取函数拟合
clear all
load('data_nasa9.mat'); load('.\molecular_w.mat') ;
%% 参数输入
c = [0.42, 0.21, 0.79, 0, 0 ,0, 0, 0, 0, 0,0 ];
m = ["H2", "O2" , "N2", "H2O" , "O" , "H", "OH", "H2O2","HO2","N","NO"];
A = [2, 0, 0, 2, 0, 1, 1, 2, 1, 0, 0; %氢原子数
     0, 2, 0, 1, 1, 0, 1, 2, 2, 0, 1; %氧原子数
     0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 1]; %氮原子数

Ma = 6;                p_1 = 101325;                
T_1 = 300;              theta_0 = deg2rad(0);     

%% 定义基本量
ns = length(m);
Mw = zeros(1,ns);  
for i = 1:ns
    Mw(i)=m_w.(m(i));
end
D = zeros(ns,9,3);     %Extract NASA 9 coefficients
for j = 1:ns
    D(j,:,:) = (coeff_nasa9.(m(j)))';
end

Mm_1 = dot(c,Mw)/sum(c); 
Rm_1 = 8.314/Mm_1;
rho_1 = p_1/(Rm_1*T_1);
a_s_1 = a_s (T_1,Mm_1,Rm_1,c,D);
u_1 = Ma * a_s_1; 

%%  激波求解
max = 200;
beta_vector_weak = zeros(1,max);  %弱解数组
theta_vector = zeros(1,max);
delta = 1;   %步进长度,对于弱解，初始时刻不低于1
for k = 2:max
    theta = theta_0 + deg2rad((k-1)*delta);
    a_0 = shock_o(u_1, p_1 ,T_1, theta, Mm_1, Mm_1, c, c, D);
    b_0 = post_reaction(a_0(1),a_0(2),c, D, A);
    e =10^(-7);
     for i = 1:max
         Mm_2 = dot(b_0,Mw)/sum(b_0); 
         a_1 = shock_o(u_1, p_1 ,T_1, theta, Mm_1, Mm_2, c, b_0, D);
         b_1 = post_reaction(a_1(1), a_1(2),c, D, A);
        if abs(b_1(1)-b_0(1))< e
            break
        end
        b_0 = b_1 + abs(b_1-b_0)*0.7;
     end
    if a_1(3)<deg2rad(10)|| a_1(3)> deg2rad(91)
        k = k-1;  
        break
    end
    beta_vector_weak(k) = rad2deg(a_1(3));
end
beta_vector_strong = zeros(1,max);  %强解数组
for p = 2:max
    theta = theta_0 + deg2rad((p-1)*delta);
    a_0 = shock_o(u_1, p_1 ,T_1, theta, Mm_1, Mm_1, c, c, D);
    b_0 = post_reaction(a_0(1),a_0(2),c, D, A);
    e =10^(-7);
     for i = 1:max
         Mm_2 = dot(b_0,Mw)/sum(b_0); 
         a_1 = shock_o_strong_solution(u_1, p_1 ,T_1, theta, Mm_1, Mm_2, c, b_0, D);
         b_1 = post_reaction(a_1(1), a_1(2),c, D, A);
        if abs(b_1(1)-b_0(1))< e
            break
        end
        b_0 = b_1 + abs(b_1-b_0)*0.7;
     end
    if a_1(3)<deg2rad(45) || a_1(3)> deg2rad(91)
        p = p-1; %#ok<*FXSET>
        break
    end
    beta_vector_strong(p) = rad2deg(a_1(3));
    theta_vector(p) = rad2deg(theta);
end
%% 绘图
k = min(k,p);
theta_vector(1)=0; beta_vector_weak(1) = 90; beta_vector_strong(1) = 90;
theta_vector(k+1)=theta_vector(k)+delta/2; 
beta_vector_weak(k+1) = (beta_vector_strong(k)+beta_vector_weak(k))/2;
beta_vector_strong(k+1)=(beta_vector_strong(k)+beta_vector_weak(k))/2 ;

plot (theta_vector(1:k+1),smooth(beta_vector_weak(1:k+1)),'b')
hold on
plot (theta_vector(1:k+1),smooth(beta_vector_strong(1:k+1)),'b')


