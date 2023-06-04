function [ U ] = W2U( W,gamma )
% tranWform W to dU
%   W: [rho u p]
%   U: [rho rho*u E]

rho = W(:,1);
u   = W(:,2);
p   = W(:,3);

U  = [rho, rho.*u, p./(gamma-1)+rho.*u.^2 /2];

end
