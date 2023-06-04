function [ W ] = U2W( U,gamma )


w   = U(:,1);
t   = U(:,2);
m   = U(:,3);

rho = w;
u   = t./w;
p   = (gamma-1).*(m-0.5.*rho.*u.^2);


W   = [rho,u,p];

end