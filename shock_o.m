%б�����ļ��㣬����
function shock_o = shock_o(u_1, p_1,T_1, theta, Mm_1, Mm_2, c ,y ,D)

Ru =8.314;
Rm_1 = Ru/Mm_1;
Rm_2 = Ru/Mm_2;
a_1 = a_s(T_1,Mm_1,Rm_1,c,D);
miu = asin(a_1/u_1);
beta_0 = theta+deg2rad(3);%��ʼ�� ����
if theta >miu
    beta_0 = theta+deg2rad(0.5);
end
max_step = 100;
e=10^(-6);%��Χ
%% ����

for i=1:max_step
    T_0 = 1 / Rm_2 *(Rm_1 *T_1 /(u_1*sin(beta_0)) +u_1*sin(beta_0)-...
    (u_1*cos(beta_0)*tan(beta_0-theta)))*(u_1*cos(beta_0)*tan(beta_0-theta));
    if T_0 < 200
        T_0 = 201;
    end
    if T_0 > 20000
        T_0 = 20000;
    end
    f = h0_m(T_0,y,Mm_2,D) +(u_1*cos(beta_0)*tan(beta_0-theta))^2 /2 - ...
        h0_m(T_1,c,Mm_1,D) - (u_1*sin(beta_0))^2 /2;
    g = -u_1^2*sin(beta_0) * cos(beta_0)+ ...
        (u_1*cos(beta_0)*tan(beta_0-theta))* u_1*(cos(beta_0)/cos(beta_0-theta)^2 - sin(beta_0)* tan(beta_0-theta))+...
        Cp_m(T_0,y,Mm_2,D) *...
        ((1/Rm_2* (u_1*cos(beta_0)*tan(beta_0-theta))*(u_1*cos(beta_0)-Rm_1*T_1/u_1*cos(beta_0)/sin(beta_0)^2))+...
        1/Rm_2*(Rm_1*T_1/(u_1*sin(beta_0))+ u_1*sin(beta_0) -2 * (u_1*cos(beta_0)*tan(beta_0-theta)))*...
        (u_1*(cos(beta_0)/cos(beta_0-theta)^2 - sin(beta_0)* tan(beta_0-theta))));
        beta_1=beta_0-f/g;
        
    if abs(beta_1-beta_0)<e
        break
    end
    beta_0=beta_1;
    
end
%% ���
rho_1 = p_1/(Rm_1*T_1);
T_2 = T_0;
u_2 = (u_1*cos(beta_0)*tan(beta_0-theta))/sin(beta_0-theta);
rho_2 = rho_1 *(u_1*(sin(beta_0)))/(u_2*sin(beta_0-theta));
p_2 = rho_2 * Rm_2 * T_2;
shock_o = [T_2, p_2, beta_0];
