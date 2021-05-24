
function t_vec = GC_solve(M,P,Q,pw_m,pz)

%% Solve the Green-Christoffel equation to get the polarization vector and 
%% the slowness vector
%% Before call this function, makesure that pz is the eigenvalue of GC equation
%% Otherwise this subroutine won't give correct result
%% Input: M,P,Q,pw_m four matrices parts of the polynomial form of GC, 
%%        M is the coefficient matrix of pz^2, P is the coefficient matrix 
%%        of pz, Q-pw_m is the constant matrix
%% Input: pz is the eigenvalue of GC equation
%% Output: t_vec: polarization vector


t_vec = zeros(6,3);


p3 = pz(1);
GC = p3^2*M + p3*P + Q - pw_m;
GC = rand_othgtra(GC);
[~,R] = qr(GC); % take QR domposition of G-C equation
A = R(1:2,2:3); % see notes
b = -[R(1,1);0];
x = A\b;
t = [1;x];
t_vec(1,:) = t/norm(t);

for i = 2:numel(pz)

    if  abs(pz(i) - pz(i-1)) > 1e-4
        p3 = pz(i);
        GC = p3^2*M + p3*P + Q - pw_m;
        GC = rand_othgtra(GC);
        [~,R] = qr(GC); % take QR domposition of G-C equation
        A = R(1:2,2:3); % see notes
        b = -[R(1,1);0];
        x = A\b;
        t = [1;x];
        t_vec(i,:) = t/norm(t);
        
    elseif abs(pz(i) - pz(i-1)) < 1e-4 % Degenerate state
        p3 = pz(i);
        GC = p3^2*M + p3*P + Q - pw_m;
        GC = rand_othgtra(GC);
        [~,R] = qr(GC);
        t = [2;1;(-2*R(1,1)-R(1,2))/R(1,3)]; % Solve the first eigenvector
        t_vec(i,:) = t/norm(t);
        
        t2 = [1;2;(-R(1,1)-2*R(1,2))/R(1,3)]; % Solve another eigenvector
        t2 = t2 - dot(t2,t)/dot(t,t)*t; % Gram-Schmidt Orthogonalization
        t_vec(i-1,:) = t2/norm(t2);
      
    elseif abs(pz(i) - pz(i+1)) < 1e-4
        continue;
    end
    
end