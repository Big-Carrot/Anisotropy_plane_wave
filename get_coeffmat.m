function [B0,B1,B2] = get_coeffmat(c, p1, p2)
%% Give the stress tensor cijkl*kj*kl  
%% px,py are known, only pz is unknow. Divide the matrix into three parts, B2*pz^2 + B1*pz + B0
%% Input: c->cijkl, elastic constant, p1,p2, px,py. rho density
%% Output: B0 const for the polynomial, B1, linear coefficient, B2 quadrant coefficient 
B0 = zeros(3,3);
B1 = zeros(3,3);
B2 = zeros(3,3);

for i = 1:3
    for j = 1:3
   B0(i,j) = c(i,1,j,1)*p1^2+c(i,1,j,2)*p1*p2+c(i,2,j,1)*p2*p1+c(i,2,j,2)*p2^2;
   B1(i,j) = c(i,1,j,3)*p1+c(i,2,j,3)*p2+c(i,3,j,2)*p2+c(i,3,j,1)*p1;
   B2(i,j) = c(i,3,j,3);
    end
end


end