%化学平衡计算
%误差控制数组、收敛因子的更改需要谨慎，容银引起较大偏差
function post_reaction = post_reaction(T, p, c, D, A)
Ru = 8.314;
p = p/101325;
a_1 = A(1,:);  a_2 = A(2,:);  a_3 = A(3,:);
max = 100;
e=[10^(-6), 10^(-8), 10^(-10), 10^(-7)];%误差控制数组
%指定迭代初值，假设已经开始反应了。
 if (c(2)-c(1)/4) >= 0.001
     y = [c(1)/2, c(2)-c(1)/4, c(3), c(4)+c(1)/2+e(4), c(5)+e(4), c(6)+e(4), c(7)+e(4), c(8)+e(4), c(9)+e(4), c(10)+e(4),c(11)+e(4)];
 else
     y = [c(1), c(2), c(3), c(4), c(5)+e(4),  c(6)+e(4), c(7)+e(4), c(8)+e(4), c(9)+e(4), c(10)+e(4),c(11)+e(4)];  
 end
 
for i = 1:max
%% 方程求解
     d = y.*(G0_i(T,D)/(Ru * T) + log(y*p/sum(y)));
     E = [dot(a_1,y), dot(a_2, y), dot(a_3, y), 0;...
         sum(a_1.* y .* a_1), sum(a_2.* y .* a_1), sum(a_3.* y .* a_1), dot(a_1,y);...
         sum(a_1.* y .* a_2), sum(a_2.* y .* a_2), sum(a_3.* y .* a_2), dot(a_2,y);...
         sum(a_1.* y .* a_3), sum(a_2.* y .* a_3), sum(a_3.* y .* a_3), dot(a_3,y)];
     %co = cond(A);
     F = [sum(d); sum(a_1.*(c+d));sum(a_2.*(c+d));sum(a_3.*(c+d))];
     L = gmres(E,F,[],e(2));  %L = [l_1, l_2, l,3, sum(x)/sum(y)] l_1, l_2是对应氢、氧方程的拉格朗日乘子
     %gmres算法一定要指定初值
     x = zeros(1,length(c));
     for j = 1:length(c)
         x(j) = -d(j) + y(j) * (L(4)+ dot(L(1:3,:),A(:,j)));
     end
%       x = real(x);
 %% 收敛因子
     omega = zeros(1,length(c));
     for j = 1:length(c)
         if x(j)< -e(1) 
             omega(j) = (e(2)-y(j))/(x(j)-y(j));
         elseif x(j) > 0
             omega(j) = 1;
         else 
             omega(j) = 1;
             x(j) = e(3);
         end
     end
     omega_k = min(omega);
     for k =1:max
         z = y + omega_k *(x - y);
         g = z.* (G0_i(T,D)/(Ru * T) +log(z*p/sum(z)));

         if sum(g) < sum(d)
             break
         end
         omega_k = omega_k*0.99;
     end
%% 方程求解
     if abs(z-y)<e(1)
         break
     end
     y = real(z);

end
post_reaction = y;