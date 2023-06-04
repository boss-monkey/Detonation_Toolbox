%% Ò»Î¬ÎÞÕ³·´Ó¦Á÷£¨±¬ºä²¨¹Ü£©¼ÆËã£¬»ùÓÚÊ±¼äÏî·ÖÁÑµÄ½âñîËã·¨
%author :Boss Monkey
%Email : Baiht0201@nuaa.edu.cn 
%²Î¿¼ÎÄÏ× Á¿ÈÈÍêÈ«ÆøÌå/»¯Ñ§·ÇÆ½ºâÆøÌåÁ÷¶¯ ALEÓÐÏÞÌå»ý¼ÆËã·½·¨ÑÐ¾¿
         %A Quasi-Steady-State Solver for the Stiff Ordinary Differential Equations of Reaction Kinetics
         %A Solver for the Stiff Ordinary Differential Equations of Reaction Kinetics
%%
clear all %#ok<CLALL>
load('data_nasa9.mat')
load('.\molecular_w.mat') ;
%% Import data
[reaction,mix] = data_import_12_23(m_w,coeff_nasa9);
m = mix.m; Mw =mix.Mw;D= mix.D;c_0 = mix.c_0;R_i = mix.R_i; ns= mix.ns;
Ru = 8.314;

%% Initializate
p_0 = 0.066*101325;  p_1 = 0.362*101325;
T_0 = 298.0;   T_1 = 621;

y_0= c_0 /dot(c_0,Mw).*Mw ;
y_1 = y_0;
rho_0 = p_0/(sum(R_i.*y_0) * T_0);
rho_1 = p_1/(sum(R_i.*y_1) * T_1);

gamma_0 = cal_gamma_T(T_0,y_0,mix);
gamma_1 = cal_gamma_T(T_1,y_1,mix);

%% Flow settings
CFL = 0.5;
time = 220*1e-6;     %total time
a = 0; b = 0.25;
dx =  0.25*1e-3;          
N = round((b-a)/dx)+1;    %number of nodes
x = linspace(a,b,N)';
x0 = 0.01;
x1 = 0.25;

%% Initializate the rho u p gamma;
rho0 =  rho_0.*(x<x0) + rho_1.*(x>=x0 & x<=x1)+rho_0.*(x>x1);    
u0 =   0.*(x<x0) + -476 *(x>=x0 & x<=x1)+0.*(x>x1);
p0 =  p_0.*(x<x0) + p_1 .*(x>=x0 & x<=x1)+p_0.*(x>x1);
gamma0 =  gamma_0.*(x<x0) + gamma_1.*(x>=x0 & x<=x1)+gamma_0.*(x>x1);
y_c = y_0.*(x<x0) + y_1.*(x>=x0 & x<=x1)+y_0.*(x>x1);

W_0 = [rho0,u0 ,p0];
U_0 = W2U(W_0,gamma0);

%% Calculation settings
current_time = 0;
maxSteps     = 1e4;
flag = 1;
step = 0;
U = U_0;
Qc = zeros(N,ns+3);

W = U2W(U,gamma0);
%% time stepping
tic
while(flag)
    
    dt = CFL*dx / max(abs(W(:,2))+sqrt(gamma0.*W(:,3)./W(:,1)));
    step = step+1;
    
    if step > maxSteps          % stoping criteria
        disp('maxSteps reached');
        break;
    end
    
    if current_time+dt >= time   % stoping criteria
        dt = time-current_time;
        flag = 0;
    end
    current_time = current_time+dt  %#ok<NOPTS>
 %%%¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª Reaction step 1¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª%%%    
    for i = 1:N    
        Qc(i,:) = cal_Qc(dt/2, W(i,1), W(i,3), y_c(i,:), reaction,mix);
    end

    W_c = W;
    W_c(:,3)=Qc(:,1);
    gamma1 =Qc(:,2);
    y_c = Qc(:,4:end);
    U_c = W2U(W_c,gamma1);
 %%%¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ªR-K step 1¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª%%%
    U_1 = U_c;
    F_U = Numerical_flux_weno(U_c,gamma1,N);
    U_1(3:N-2,:) = U_c(3:N-2,:) - dt/dx * F_U;
 %%%¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ªR-K step 2¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª%%%    
    U_2 = U_c;
    F_U = Numerical_flux_weno(U_1,gamma1,N);
    U_2(3:N-2,:) = 3/4*U_c(3:N-2,:)+1/4*(U_1(3:N-2,:) - dt/dx * F_U);
 %%%¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ªR-K step 3¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª%%%    
    F_U = Numerical_flux_weno(U_2,gamma1,N);
    U(3:N-2,:) = 1/3*U_c(3:N-2,:)+2/3*(U_2(3:N-2,:) - dt/dx * F_U);
 %%%¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ªBoundary value¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª%%%     
    U(1:2,:) = [U(3,:);U(3,:)];
    U(N-1:N,:) = [U(N-2,:);U(N-2,:)];
 %%%¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª Flash¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª%%%    
    gamma1 = cal_gamma_E(U,gamma1,y_c,mix);
    W_f = U2W(U,gamma1);
 %%%¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª Reaction step 2¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª%%%       
    for i = 1:N    
        Qc(i,:) = cal_Qc(dt/2, W_f(i,1), W_f(i,3), y_c(i,:), reaction,mix);
    end

    W_f (:,3)=Qc(:,1);
 %%%¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª Flash¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª¡ª%%%    
    W = W_f;
    gamma0 =Qc(:,2);
    y_c = Qc(:,4:end);
    W(1:3,2) = [0;0;0];   %Boundary conditions
end
toc
%% Result
rho = W(:,1);
u =  W(:,2);
p =  W(:,3);
T =  W(:,3)./(W(:,1).*sum(y_c.*R_i,2));

figure(1);
plot(x,p,'k-')
figure(2);
plot(x,y_c(:,1),'k-')
figure(3);
plot(x,T,'k-')
figure(4);
plot(x,u,'k-')

y_rho = linspace(min(rho),max(rho),length(rho));
y_u = linspace(min(u),max(u),length(u));
y_p = linspace(min(p),max(p),length(p));
y_T = linspace(min(T),max(T),length(T));

z_rho = meshgrid(rho,x);
z_u = meshgrid(u,x);
z_p =  meshgrid(p,x);
z_T =  meshgrid(T,x);

figure(5);    %contour
%rho
subplot(4,1,1)
contourf(x,y_rho,z_rho,100,'linestyle','none');
xlabel('x');
ylabel('rho');
title(['dx=',num2str(dx),'  CFL=',num2str(CFL),' t=',num2str(time)]);
% u
subplot(4,1,2)
contourf(x,y_u,z_u,100,'linestyle','none');
xlabel('x');
ylabel('u');

% pressure
subplot(4,1,3)
contourf(x,y_p,z_p,100,'linestyle','none');
xlabel('x');
ylabel('p');

subplot(4,1,4)
contourf(x,y_T,z_T,100,'linestyle','none');
xlabel('x');
ylabel('T');
colormap(jet);