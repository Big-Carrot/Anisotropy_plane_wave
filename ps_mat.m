function Ps = ps_mat(omega, p_vec, h)
%% Propagation matrix
%% Input:p_vec, the slowness of all waves. 6*3 matrix, 6 for 3 upgoing 
%%       3 downgoing waves, 3 for xyz components
%% Input:h, layer thickness
%% Output: Ps, phase-shifting matrix in this layer

Ps = eye(6);
Ps(1,1) = exp(1i*p_vec(1,3)*h*omega);
Ps(2,2) = exp(1i*p_vec(2,3)*h*omega);
Ps(3,3) = exp(1i*p_vec(3,3)*h*omega);
Ps(4,4) = exp(-1i*p_vec(4,3)*h*omega);
Ps(5,5) = exp(-1i*p_vec(5,3)*h*omega);
Ps(6,6) = exp(-1i*p_vec(6,3)*h*omega);


end