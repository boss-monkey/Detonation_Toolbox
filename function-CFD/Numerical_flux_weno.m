%有限差分方法WENO格式数值通量的求解
function F_U = Numerical_flux_weno(U_1,gamma,N)

e = 1e-6;
H_p = zeros(N,3);
H_n = zeros(N,3);
H_numerical_p= zeros(N,3);
H_numerical_n= zeros(N,3);
gamma_k = [1/10,3/5,3/10]';

for j = 1:N
    A1 = U2A_sgn(U_1, j , 1, gamma(j));                %A+   
    A2 = U2A_sgn(U_1, j ,-1, gamma(j));                %A-
    H_p(j,:) = A1*U_1(j,:)';
    H_n(j,:) = A2*U_1(j,:)';
end
for i = 3:N-2
    H_p_numerical = 1/6*[(2*H_p(i-2,:)-7*H_p(i-1,:)+11*H_p(i  ,:));
                         ( -H_p(i-1,:)+5*H_p(i  ,:)+ 2*H_p(i+1,:));
                         (2*H_p(i  ,:)+5*H_p(i+1,:)-   H_p(i+2,:))];

    H_n_numerical = 1/6*[(2*H_n(i+2,:)-7*H_n(i+1,:)+11*H_n(i  ,:));
                          (-H_n(i+1,:)+5*H_n(i  ,:)+ 2*H_n(i-1,:));
                         (2*H_n(i  ,:)+5*H_n(i-1,:)-   H_n(i-2,:))];  

    beta_kp =[(13/12*(H_p(i-2,:)-2*H_p(i-1,:)+H_p(i  ,:)).^2 + 1/4*(H_p(i-2,:)-4*H_p(i-1,:)+3*H_p(i ,:)).^2);
              (13/12*(H_p(i-1,:)-2*H_p(i  ,:)+H_p(i+1,:)).^2 + 1/4*(H_p(i-1,:)-H_p(i+1,:)).^2);
              (13/12*(H_p(i  ,:)-2*H_p(i+1,:)+H_p(i+2,:)).^2 + 1/4*(H_p(i+2,:)-4*H_p(i+1,:)+3*H_p(i ,:)).^2)];
    w_kp_ = gamma_k./(e+beta_kp).^2;
    w_kp = w_kp_./sum(w_kp_);

    beta_kn =[(13/12*(H_n(i+2,:)-2*H_n(i+1,:)+H_n(i  ,:)).^2 + 1/4*(H_n(i+2,:)-4*H_n(i+1,:)+3*H_n(i ,:)).^2);
              (13/12*(H_n(i-1,:)-2*H_n(i  ,:)+H_n(i+1,:)).^2 + 1/4*(H_n(i-1,:)-H_n(i+1,:)).^2);
              (13/12*(H_n(i  ,:)-2*H_n(i-1,:)+H_n(i-2,:)).^2 + 1/4*(H_n(i-2,:)-4*H_n(i-1,:)+3*H_n(i ,:)).^2)];
    w_kn_ = gamma_k./(e+beta_kn).^2;
    w_kn = w_kn_./sum(w_kn_);

    H_numerical_p(i,:) = sum(H_p_numerical.* w_kp);
    H_numerical_n(i,:) = sum(H_n_numerical.* w_kn);


end

H_numerical_p(1:2,:)   = [H_numerical_p(3,:);H_numerical_p(3,:)];
H_numerical_p(N-1:N,:) = [H_numerical_p(N-2,:);H_numerical_p(N-2,:)];
H_numerical_n(1:2,:)   = [H_numerical_n(3,:);H_numerical_n(3,:)];
H_numerical_n(N-1:N,:) = [H_numerical_n(N-2,:);H_numerical_n(N-2,:)];

F_U = H_numerical_p(3:N-2,:)-H_numerical_p(2:N-3,:)+...    %F_p
      H_numerical_n(4:N-1,:)-H_numerical_n(3:N-2,:);       %F_n
    
end
