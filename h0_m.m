%混合气单位质量焓
function h0_m = h0_m(T, c, Mm ,D)
data_cal = data_cal_vector(T,D);
h0 =  data_cal(2,:);
h0_m = dot(c,h0)/sum(c)/Mm;
end
