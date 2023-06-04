function F_U = Numerical_flux_upwind(U_1,gamma,N)
e = 1e-16;
F_U = zeros(N-4,3);
%TVD limiters  
r_r1 = (U_1(3:N-2,:)-U_1(2:N-3,:))./ (U_1(4:N-1,:)-U_1(3:N-2,:)+e);
phi_r1 = 2*r_r1./(1+r_r1.^2);
r_r2 = (U_1(2:N-3,:)-U_1(1:N-4,:))./ (U_1(3:N-2,:)-U_1(2:N-3,:)+e);
phi_r2 = 2*r_r2./(1+r_r2.^2);
r_l1 = (U_1(4:N-1,:)-U_1(5:N  ,:))./ (U_1(3:N-2,:)-U_1(4:N-1,:)+e);
phi_l1 = 2*r_l1./(1+r_l1.^2);
r_l2 = (U_1(3:N-2,:)-U_1(4:N-1,:))./ (U_1(2:N-3,:)-U_1(3:N-2,:)+e);
phi_l2 = 2*r_l2./(1+r_l2.^2); 
%R-K 1step
for j = 3:N-2
    Al = U2A_sgn(U_1, j ,-1, gamma(j));                %A-
    Ar = U2A_sgn(U_1, j , 1, gamma(j));                %A+
    F_U(j-2,:) = (Al*(U_1(j+1,:)+phi_l1(j-2,:).*(U_1(j  ,:)-U_1(j+1  ,:))/2 - ...
                      U_1(j  ,:)-phi_l2(j-2,:).*(U_1(j-1,:)-U_1(j  ,:))/2 )'...
                 +Ar*(U_1(j  ,:)+phi_r1(j-2,:).*(U_1(j+1,:)-U_1(j  ,:))/2 - ...
                      U_1(j-1,:)-phi_r2(j-2,:).*(U_1(j  ,:)-U_1(j-1,:))/2 )')';
end
end
