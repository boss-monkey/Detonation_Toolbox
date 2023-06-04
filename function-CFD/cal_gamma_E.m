%通过牛顿方法根据气体内能（守恒量）求解比热比
function gamma = cal_gamma_E(U,gamma0,y_i,mix)

    u =  U(:,2)./ U(:,1);
    e =  U(:,3)./ U(:,1) - 0.5 * u.^2;
    T0 = e.*(gamma0-1)./sum(y_i.*mix.R_i,2);
    N =length(T0);
    data_cal_st = data_cal_vector(298.15,mix.D);                  %标准温度下热力学系数
    h_i0 = (data_cal_st(2,:) - data_cal_st(1,:)*298.15)./mix.Mw;   %生成焓
    h_i = zeros(N,mix.ns);
    cv_i = zeros(N,mix.ns);
    
for j = 1:100    
    for i = 1:N  
        data_cal = data_cal_vector(T0(i),mix.D);
        h_i(i,:) = data_cal(2,:)./mix.Mw;
        cv_i(i,:) = data_cal(1,:)./mix.Mw - mix.R_i;
    end
    h_star = h_i - h_i0 ;                                         % 热焓
    f = sum(y_i.*h_star,2)-sum(y_i.*mix.R_i,2).*T0 -e;
    g = sum(y_i.*(cv_i),2);
    T = T0 -f./g;
    if max(abs(T-T0))<1e-6
        break
    end
    T0 = T; 
    
end
   gamma = (e+T.*sum(y_i.*mix.R_i,2))./e;
  

end
