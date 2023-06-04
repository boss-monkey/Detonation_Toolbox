%混合气单位质量比热
function Cp_m = Cp_m(T, c, Mm ,D)
data_cal = data_cal_vector(T,D);
cp = data_cal(1,:);
Cp_m = dot(cp,c)/sum(c)/Mm;
end
