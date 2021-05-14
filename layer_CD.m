function [C,D] = layer_CD(k_vec,t_vec, Cij)
%% Calculate the composition matrix [Up,Dw]->[u,t] 
%% and decomposition matrix [u,t]->[Up,Dw] 
%% Input: k_vec, 6*3 matrix, wavenumber for all waves.3 for x,y,z components, 6 for 3 upgoing,
%%        3 downgoing waves.
%% Input: t_vec, 6*3 matrix, polarization vector for all waves. Similar to k_vec
%% Input: omega: Angular frequency w = 2*pi*f
%% Input: Cij, Cij is the elastic tensor.
%% Output: C matrix, composition matrix, C = [U;W*Cij*E](See notes)
%% Output: D matrix, the inverse of C matrix, decompisition matrix
%%         D = inv(C)

k = k_vec;
t = t_vec;
%% Calculate displacement matrix U
U = transpose(t);

%% Calculate the C matrix

 % E matrix
 E = zeros(6,6);
 E(1,:) = 1i* t(:,1).*k(:,1);
 E(2,:) = 1i* t(:,2).*k(:,2);
 E(3,:) = 1i* t(:,3).*k(:,3);
 E(4,:) = 1i*(t(:,3).*k(:,2)+t(:,2).*k(:,3));
 E(5,:) = 1i*(t(:,3).*k(:,1)+t(:,1).*k(:,3));
 E(6,:) = 1i*(t(:,1).*k(:,2)+t(:,2).*k(:,1));
 
 % w matrix
 w = zeros(3,6);
 w = [0 0 0 0 1 0;
      0 0 0 1 0 0;
      0 0 1 0 0 0];
 
 Dp = w*Cij*E;
 Up = U;
 
 C = [Up;Dp];
 
 D = pinv(C);
 
end
 



end