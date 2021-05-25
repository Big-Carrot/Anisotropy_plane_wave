function t_vec = GC_solveQr(M,P,Q,pw_m,pz)

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

[pz_u,~,ic] = unique(real(pz),'legacy');

for i = 1:numel(pz_u)
    
    ind = find(ic == i);
    p3 = pz_u(i);
    
    
    if numel(ind) == 1 %non-degenerate case
        GC = p3^2*M + p3*P + Q - pw_m;
        [GC,O] = rand_othgtra(GC);
        [~,R] = qr(GC); % take QR domposition of G-C equation
        A = R(1:2,2:3); % see notes
        b = -[R(1,1);0];
        x = A\b;
        t = [1;x];
        t_vec(ind,:) = O*t/norm(t);
        
    else   % degenerate case
        GC = p3^2*M + p3*P + Q - pw_m; 
        [GC,O] = rand_othgtra(GC);
        [~,R] = qr(GC); % take QR domposition of G-C equation
        t = [2;1;(-2*R(1,1)-R(1,2))/R(1,3)]; % Solve the first eigenvector
        t1 = t/norm(t);
        
        t = [1;2;(-R(1,1)-2*R(1,2))/R(1,3)]; % Solve another eigenvector
        t = t - dot(t,t1)/dot(t1,t1)*t1; % Gram-Schmidt Orthogonalization
        t2 = t/norm(t);
        
        t1 = O*t1;
        t2 = O*t2;
        
        t_vec(ind,:) = [t1';t2'];
        
    end
    
    
end




end