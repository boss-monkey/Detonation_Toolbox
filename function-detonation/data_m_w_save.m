%% H_2 page130
m_w.H2 = 2.01588*10^(-3);
%% 0_2 page166
m_w.O2= 31.99880*10^(-3);
%% N_2 page156
m_w.N2= 28.01340*10^(-3);
%% H page128
m_w.H= 1.00794*10^(-3);
%% HNO page129
m_w.HNO= 31.01404*10^(-3);
%% HNO_2 page129
m_w.HNO2= 47.01344*10^(-3);
%% HO_2 page130
m_w.HO2= 33.00674*10^(-3);
%% H_2O page131
m_w.H2O= 18.01528*10^(-3);
%% H_2O_2 page131
m_w.H2O2= 34.01468*10^(-3);
%% N page153
m_w.N= 14.00670*10^(-3);
%% NH page154
m_w.NH= 15.01464*10^(-3);
%% NH_2 page154
m_w.NH2= 16.02258*10^(-3);
%% NH3 page154
m_w.NH3= 17.03052*10^(-3);
%% NO page155
m_w.NO= 30.00610*10^(-3);
%% NO_2 page155
m_w.NO2= 46.00550*10^(-3);
%% O page164
m_w.O= 15.99940*10^(-3);
%% OH page165
m_w.OH= 17.00734*10^(-3);
%% O_3 page166
m_w.O3= 49.99820*10^(-3);
%% CO page166
m_w.CO= 28.0101000*10^(-3);
%% CO page166
m_w.CO2= 44.0095000*10^(-3);
%% Ar page55
m_w.Ar = 39.9480000*10^(-3);
%% save
%M_w = struct()
save molecular_w.mat m_w
load('.\molecular_w.mat') ;