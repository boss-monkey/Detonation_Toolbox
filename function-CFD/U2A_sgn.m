function [A_sgn] = U2A_sgn(U, j, sgn, gamma)
% U: [rho, rho*u, E]
% global gamma    %È«¾Ö±äÁ¿½µµÍËÙ¶È

rho = U(j,1);
u   = U(j,2)/U(j,1);
p   = (gamma-1)*(U(j,3)-0.5.*rho.*u.^2);

a =  (gamma*p/rho)^0.5;
H =  (U(j,3)+p)/rho;


K = [  1        1      1 
      u-a       u     u+a
     H-u*a     0.5*u^2   H+u*a];
 
L = [   u-a        0       0
        0        u     0
        0        0      u+a];       % Eigenvalues

 
A_sgn = K*((L+sgn.*(L.^2+1e-12).^0.5)/2)/K;

end
