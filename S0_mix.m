%µ¥Î»Ä¦¶ûìØ
function S0_mix = S0_mix(T, c, D)
data_cal = data_cal_vector(T,D);
s0 = data_cal(3,:);
S0_mix = dot(s0,c)/sum(c);
end

