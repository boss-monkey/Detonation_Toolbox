%% 爆轰波ZND结构求解，基于ODE15
%author :Boss Monkey
%Email : Baiht0201@nuaa.edu.cn 
%%
clear all %#ok<CLALL>
load('data_nasa9.mat')
load('.\molecular_w.mat') ;
%% 数据导入
[reaction,mix] = data_import(m_w,coeff_nasa9);
m = mix.m; Mw =mix.Mw;D= mix.D;c_0 = mix.c_0; ns= mix.ns;
Ru = 8.314;

%% 初值
p_0 = 0.362*101325;
T_0 = 298.15;
u_0 = 2000;%1968.112805277358;
Mm_0 = dot(c_0,Mw)/sum(c_0); 
Rm_0 = Ru/Mm_0;
rho_0 = p_0/(Rm_0 * T_0);
%% 冯纽曼条件（诱导正激波后参数）
s_1 = shock_n(u_0, p_0 ,T_0, Mm_0, Mm_0, c_0, c_0, D);    
p_1 = s_1(2);
T_1 = s_1(1);
u_1 = s_1(3);
rho_1 = p_1/(Rm_0*T_1);
%% 基本量计算
R_i = Ru./ Mw;
y_i = zeros(1,ns);   %各组分质量分数
for i = 1:ns
    y_i(i)=c_0(i)/dot(c_0,Mw) *Mw(i);
end
ma = u_1*rho_1;  %mass flux
mo = u_1^2*rho_1 + p_1;
%% 求解微分方程
flux = [ma, mo];
y0 = [u_1,y_i]';
tel = [0 0.004];
options = odeset('RelTol',1.e-5,'AbsTol',1.e-12,'Stats','on');
t0 = cputime;

out = ode15s(@ode_ZND,tel,y0,options, reaction ,mix,flux,D);
disp(['CPU time = ' num2str(cputime - t0)]);
%% 输出
yt =  out.y(2:end,:);
Rt = sum(R_i'.*yt);
rhot = ma./out.y(1,:);
pt = mo - rhot.* out.y(1,:).^2;
Tt = pt./rhot./Rt;

figure(1);
plot(out.x,Tt,'-k');
xlabel('x/m');
ylabel('Temperature');
title(['Final T = ' num2str(Tt(:,end)) ' K']);
hold on
figure(2);
plot(out.x,pt,'-k');
xlabel('x/m');
ylabel('P/pa');
hold on
figure(3);
plot(out.x,out.y(1,:),'-k');
xlabel('x/m');
ylabel('u/m s-1');
out.y(1,end)