%% 零维燃烧（等容爆炸）计算，基于alpha-QSS方法
%author :Boss Monkey
%Email : Baiht0201@nuaa.edu.cn 
%参考文献 A Quasi-Steady-State Solver for the Stiff Ordinary Differential Equations of Reaction Kinetics
         %A Solver for the Stiff Ordinary Differential Equations of Reaction Kinetics
%%
clear all %#ok<CLALL>
load('data_nasa9.mat')
load('.\molecular_w.mat');
%% Import data
[reaction,mix] = data_import(m_w,coeff_nasa9);
m = mix.m; Mw =mix.Mw;c_0 = mix.c_0;D = mix.D ;ns= mix.ns;
Ru = 8.314;
%% Initializate
p_0 = 101325;
T_0 = 1000;
T = T_0;
p = p_0;
%% Calculate the basic quantity
y_i= c_0 /dot(c_0,Mw).*Mw ;
%% Reaction
tic
t= 1e-3;
cal = cal_aQSS(t, T, p, y_i, reaction,mix);
toc
