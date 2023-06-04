function gamma  = cal_gamma_T(T,y_i,mix)
    
    data_cal_st = data_cal_vector(298.15,mix.D);                  %标准温度下热力学系数
    data_cal = data_cal_vector(T,mix.D);
    
    h_i = data_cal(2,:)./mix.Mw;
    h_i0 = (data_cal_st(2,:) - data_cal_st(1,:)*298.15)./mix.Mw;   %生成焓
    h_star = h_i - h_i0 ;                                          % 热焓
    e_star =h_star  - mix.R_i*T;
    gamma = sum(h_star .*y_i)/sum(e_star.*y_i);

end
