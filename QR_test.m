clear all
m = 6;
n = m+1;

A = rand(m,n);
A(m+1,:) = A(1,:) + 3*A(2,:);
% A(m+1,m+1) = A(m+1,m+1) + 1e-12; %perturbation term, check the
% sensitivity of the QR decomposition.

det(A)

[Q,R] = qr(A)
det(Q)
