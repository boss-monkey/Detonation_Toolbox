clear all
%% 初值
T = 2000;
p = 101325*10;
p_c = p/101325;
R_u = 8.3144;
load('data_nasa9.mat')
load('.\molecular_w.mat') ;
% H_2, O_2, H_2O, O, H, OH, H_2O_2
c_0 = [2, 1, 3.76, 0, 0 ,0, 0, 0, 0, 0, 0 ];
m = ["H2", "O2" , "N2", "H2O" , "O" , "H", "OH", "H2O2","HO2","N","NO"];
A = [2, 0, 0, 2, 0, 1, 1, 2, 1, 0, 0; %氢原子数
     0, 2, 0, 1, 1, 0, 1, 2, 2, 0, 1; %氧原子数
     0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 1]; %氮原子数
ns = length(c_0);
Mw = zeros(1,ns);  
for i = 1:ns
    Mw(i)=m_w.(m(i));
end
D = zeros(ns,9,3);     %Extract NASA 9 coefficients
for j = 1:ns
    D(j,:,:) = (coeff_nasa9.(m(j)))';
end
R_i = 8.31442./ Mw; 

x=post_reaction(T, p, c_0, D, A);
mix = struct('m',m,'Mw',Mw,'D',D,'c_0',c_0,'R_i',R_i,'ns',ns);
cal_gamma_T(T,x,mix)

x_mole_fraction = x/sum(x);


