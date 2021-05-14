function Ps = ps_mat(k_vec, h)
%% Propagation matrix
%% Input:k_vec, the wavenumber of all waves. 6*3 matrix, 6 for 3 upgoing 
%%       3 downgoing waves, 3 for xyz components
%% Input:h, layer thickness
%% Output: Ps, phase-shifting matrix in this layer

Ps = eye(6);
Ps(1,1) = exp(1i*k_vec(3,1)*h);
Ps(2,2) = exp(1i*k_vec(3,2)*h);
Ps(3,3) = exp(1i*k_vec(3,3)*h);
Ps(4,4) = exp(1i*k_vec(3,4)*h);
Ps(5,5) = exp(1i*k_vec(3,5)*h);
Ps(6,6) = exp(1i*k_vec(3,6)*h);


end