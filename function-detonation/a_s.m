function a_s = a_s (T,Mm,Rm,c,D)
gamma = Cp_m(T,c,Mm,D)/(Cp_m(T,c,Mm,D)-Rm);
a_s = (gamma* Rm * T)^0.5;
end
