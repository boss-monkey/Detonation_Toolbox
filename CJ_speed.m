%% 求C-J速度，注意收敛因子的更改需谨慎
clear all;
load('data_nasa9.mat')
load('.\molecular_w.mat') ;
%% Import
p_1 = 101325;   T_1 = 298;

c = [2, 1, 3.76, 0, 0 ,0, 0, 0, 0, 0, 0 ];
m = ["H2", "O2" , "N2", "H2O" , "O" , "H", "OH", "H2O2","HO2","N","NO"];
A = [2, 0, 0, 2, 0, 1, 1, 2, 1, 0, 0; %氢原子数
     0, 2, 0, 1, 1, 0, 1, 2, 2, 0, 1; %氧原子数
     0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 1]; %氮原子数
%% basic data
ns = length(m);
Mw = zeros(1,ns);  
for i = 1:ns
    Mw(i)=m_w.(m(i));
end
D = zeros(ns,9,3);     %Extract NASA 9 coefficients
for j = 1:ns
    D(j,:,:) = (coeff_nasa9.(m(j)))';
end
Ru = 8.314;
Mm_1 = dot(c,Mw)/sum(c); 
Rm_1 =Ru/Mm_1;
rho_1 = p_1/(Rm_1*T_1); 

%% solve
u_0 = 1500;
a_0 = shock_n(u_0, p_1 ,T_1, Mm_1, Mm_1, c, c, D);
b_0 = post_reaction(a_0(1),a_0(2),c, D, A);

e = 1e-1;
h = zeros(); 
tic
for i = 1:50
    T_min = 1000;
    T_max = 5000;
    for j = 1:50 
        Mm_2 = dot(b_0,Mw)/sum(b_0);
        f = h0_m(a_0(1),b_0,Mm_2,D)-h0_m(T_1,c,Mm_1,D);
        g = (a_0(2)-p_1)*((Ru/Mm_2)*a_0(1)/a_0(2)+Rm_1*T_1/p_1)/2;
        if abs(f-g)<e
            break
        end
        if (f-g)>e
            T_max = a_0(1);
            a_0(1)=(a_0(1)+T_min)/2;
        elseif (f-g)<-e 
            T_min = a_0(1);
            a_0(1)=(a_0(1)+T_max)/2;
        b_0 = post_reaction(a_0(1),a_0(2),c, D, A);
        end
    end

    Rm_2 = Ru/Mm_2;
    gamma_2 = Cp_m(a_0(1),b_0,Mm_2,D)/(Cp_m(a_0(1),b_0,Mm_2,D)-Rm_2);
    r = (a_0(2)-p_1)/(gamma_2*a_0(2))+1;
    p_k = p_1* r *(Rm_2/Rm_1)*(a_0(1)/T_1);
    if abs(a_0(2)- p_k)<e
        break
    end
    x = 0.4;   %收敛因子
    a_0(2)= a_0(2)*x +(1-x)*p_k;
    b_0 = post_reaction(a_0(1),a_0(2), c, D, A);
    
    h(i) = a_0(1);
end
toc
%% output
p_2 = a_0(2);
T_2 = a_0(1);
rho_2 = rho_1*r;
u_1 = (gamma_2*p_2*rho_2/rho_1^2)^0.5;
u_2 = u_1/r;
c_2 = (gamma_2*Rm_2*T_2)^0.5;
plot(h) 
hold on
%% check
q = h0_m(T_1,c,Mm_1,D)- h0_m(T_1,b_0,Mm_2,D);  %标准摩尔反应热
a1 = h0_m(T_2,b_0,Mm_2,D)- h0_m(T_1,b_0,Mm_2,D);
a2 = 1/2*(u_1^2 - u_2^2) + q;

if abs(a1-a2)<e
    fprintf('计算收敛');
else
    fprintf('计算未收敛');
end

