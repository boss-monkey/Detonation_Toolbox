%% 零维燃烧（等容爆炸）计算，基于ODE15
%author :Boss Monkey
%Email : Baiht0201@nuaa.edu.cn 
%%
clear all %#ok<CLALL>
load('data_nasa9.mat')
load('.\molecular_w.mat') ;
%% Import data
[reaction,mix] = data_import(m_w,coeff_nasa9);
m = mix.m; Mw =mix.Mw;D= mix.D;c_0 = mix.c_0; ns= mix.ns;
Ru = 8.314;
%% 初值
p_0 = 101325;
T_0 = 950;
y_i= c_0 /dot(c_0,Mw).*Mw ;
R_i = Ru./ Mw;
Rm_0 = dot(R_i,y_i);
rho = p_0/(Rm_0 * T_0);

%% 求解微分方程
tic
y0 = [T_0,y_i]';
tel = [0 0.001];
options = odeset('RelTol',1.e-5,'AbsTol',1.e-12);
t0 = cputime;

out = ode15s(@ode_reaction,tel,y0,options, reaction ,mix, rho);
disp(['CPU time = ' num2str(cputime - t0)]);

toc
%% plot the temperature and O mole fractions.
figure(1);
plot(out.x,out.y(1,:));
xlabel('time');
ylabel('Temperature');
title(['Final T = ' num2str(out.y(1,end)) ' K']);
hold on
figure(2);
plot(out.x,out.y(5,:));
xlabel('time');
ylabel('O');
out_end = out.y(:,end)';
p = rho*dot(out_end(2:end),R_i)*out_end(1);
more_farc = rho*out_end(2:end)./Mw/sum(rho*out_end(2:end)./Mw);

