function A_trans = rand_othgtra(A)
%% generate a random orthogonal transformation matrix O*Transpose(O) = I
%% Input: m*m matrix A
%% Output: transformed matrix A_trans

[m,n] = size(A);
if m~=n
    disp('Error, A should be a square matrix');
    return
else
[O,~] = qr(rand(m,m));
OT = transpose(O);
A_trans = O*A*OT;
end
end