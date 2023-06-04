function [reaction,mix] = data_import_12_23(m_w,coeff_nasa9)

[~,m] = xlsread('H2-Air_12_23.xlsx','B1:M1');  %物质名称
c_0 = xlsread('H2-Air_12_23.xlsx','B27:M27');  %物质摩尔数
v_1 = xlsread('H2-Air_12_23.xlsx','B3:M25');   %反应物系数矩阵
v_2 = xlsread('H2-Air_12_23.xlsx','O3:Z25');   %生成物系数矩阵
A = xlsread('H2-Air_12_23.xlsx','AB3:AB25');      %指前因子
beta = xlsread('H2-Air_12_23.xlsx','AC3:AC25');   %温度指数
Ea = xlsread('H2-Air_12_23.xlsx','AD3:AD25');     %活化能
Mr = xlsread('H2-Air_12_23.xlsx','AG3:AR8');
G = xlsread('H2-Air_12_23.xlsx','AF3:AF25');
m = string(m);
fprintf('read files down'); 
nr = length(A);
nt = sum(G);
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
reaction =  struct('v_1',v_1,'v_2', v_2,'A', A,'beta',beta,'Ea',Ea,'Mr', Mr, 'nr',nr,'nt' ,nt);
mix = struct('m',m,'Mw',Mw,'D',D,'c_0',c_0,'R_i',R_i,'ns',ns);
