%等熵膨胀到指定压力
function isentropic_cal = isentropic_cal_pressure(p_1, p_2, T_1 ,u_1, c, y, Mw, D)
Ru = 8.314;
Mm_1 = dot(c,Mw)/sum(c); 
Mm_2 = dot(y,Mw)/sum(y); 

max = 20;
e=10^(-8);%误差范围
%% 计算
T_2 = 1000; %迭代初值
for i=1:max
    f = S0_mix(T_2, y, D)/Ru- S0_mix(T_1, c, D)/Ru+log(p_1/p_2);
    g = Cp_m(T_2,y, 1 ,D)/(Ru * T_2); %M_m混合气体相对质量给1
    T_k=T_2 - f/g;
    if abs(T_k - T_2)<e
        break
    end
    T_2=T_k;
end
%% 输出
u_2 = (2*(h0_m(T_1,c,Mm_1,D) + u_1^2/2 - h0_m(T_2,y,Mm_2,D)))^0.5;
isentropic_cal = [T_2,p_2,u_2];

