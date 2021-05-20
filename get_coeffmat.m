function [B0,B1,B2] = get_coeffmat(c, k1, k2, rho)

B0 = zeros(3,3);
B1 = zeros(3,3);
B2 = zeros(3,3);

for i = 1:3
    for j = 1:3
   B0(i,j) = c(i,1,j,1)*k1^2+c(i,1,j,2)*k1*k2+c(i,2,j,1)*k2*k1+c(i,2,j,2)*k2^2 - rho;
   B1(i,j) = c(i,1,j,3)*k1+c(i,2,j,3)*k2+c(i,3,j,2)*k2+c(i,3,j,1)*k1;
   B2(i,j) = c(i,3,j,3);
    end
end


end