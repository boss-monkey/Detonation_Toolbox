function G0_i = G0_i(T,D)
data_cal = data_cal_vector(T,D);
h = data_cal(2,:);
s = data_cal(3,:);
G0_i = h - T * s;   
end