clear;
close all;
clc;
set(0, 'DefaultLineLineWidth', 1);

%% Set up the model
nlayer = 3;
for j =1: nlayer
    layer(j).vp = 6; %p wave velocity km/s
    layer(j).vs = layer(j).vp  /2;
    layer(j).rho = 2; % 1e3kg/m^3
    layer(j).gamma = 0; % Thomsen gamma
    layer(j).epsilon = 0; % Thomsen epsilon
    layer(j).delta = 0; % Thomsen delta
    layer(j).h = .500;  % layer thichness m
    
    %
    layer(1).vp = 5;
    layer(1).vs = 3;
    layer(1).h = 0.2;
    
    layer(2).vp = 12;
    layer(2).vs = 5;
    layer(2).h = 10;
    layer(2).gamma = 0.1; 
    
%   
    layer(3).vp = 5;
    layer(3).vs = 3;
    layer(3).h = 0.2;
    layer(3).rho = layer(1).rho; 
    %
    if j==1
        layer(j).z1 = 0; 
        layer(j).z2 = layer(j).h; 
    else
        layer(j).z1 = layer(j-1).z2; 
        layer(j).z2 = layer(j).z1 + layer(j).h; 
    end
    
    % set rank-4 tensor c; VTI; then rotate 
    [Cij] = thomsen2Cij_CWP(layer(j).vp, layer(j).vs, layer(j).epsilon, ...
        layer(j).delta, layer(j).gamma,layer(j).rho);
    cijkl = C2toc4(Cij);
    corig = cijkl;
    
    % rot matrix for cijkl
    theta = 40 * pi/180;
    Ry= roty(theta);
    Rz= rotz(theta);
    R = eye(3);
    if j==2
        R = Ry*Rz;
    end
    % rotate cijkl to cr_ijkl using the rotation matrix R
    cr = rot1234(corig, R);
    Cij = c4toC2(cr);
    % set layer cijkl 
    layer(j).c4 = cr;
    layer(j).C2 = Cij;  
end

%% Ray parameter set-up
% inc_ind: The incident wave polarization, 1-p, 2,3-s.
% A = [s_up,s_up,p_up,p_down,s_down,s_down]
inc_azi = 0; % degree
inc_dip = 0; % degree
A_vec = zeros(6,1);
inc_ind = 1; % P wave incidence
if inc_ind == 1
A_vec(inc_ind+3) = 1;    
p = 1/layer(1).vp;
p1 = p*sind(inc_dip)*cosd(inc_azi);
p2 = p*sind(inc_dip)*sind(inc_azi);
elseif inc_ind == 2|| inc_ind ==3 % S wave incidence
A_vec(inc_ind+3) = 1 ;    
p = 1/layer(1).vs;
p1 = p*sind(inc_dip)*cosd(inc_azi);
p2 = p*sind(inc_dip)*sind(inc_azi);    
end

%% Layer iteration for propagation matrix D,C,Ps
for j = 1:nlayer
    
 rho = layer(j).rho;
 Cij = layer(j).C2;
 px = p1;
 py = p2;
 [p_vec, t_vec] = christoffel(rho, Cij, px, py);   
 [C,D] = layer_CD(p_vec,t_vec, Cij);
 layer(j).Cmat = C;
 layer(j).Dmat = D;
 layer(j).pv = p_vec;
 layer(j).tv = t_vec;
 
end

%% Frequency iteration, get propagation matrix 
dt = 1e-3; % 1ms sampling interval
Fs = 1/dt; % Sampling frequency
T = 2; % Total time;
nfft = round(T/dt);
nfft = nfft + mod(nfft,2);
df = Fs/nfft;
ff = [0:nfft/2,-nfft/2+1:-1]*df;

ampt = zeros(6,nfft);
ampr = zeros(6,nfft);


A_dw = A_vec(4:6);
for ifre = 1:nfft/2+1
    
Pmat = eye(6,6);


omega = df*ifre*2*pi;
C1 = layer(1).Cmat;

for j = 2:nlayer-1
Cj = layer(j).Cmat; 
Dj = layer(j).Dmat;
p_vec = layer(j).pv;
h = layer(j).h;
Ps = ps_mat(omega, p_vec, h);
Pmat = Cj*Ps*Dj*Pmat;
end

Dn = layer(nlayer).Dmat;
Pmat = Dn*Pmat;


P11 = Pmat(1:3,1:3);
P12 = Pmat(1:3,4:6);
P21 = Pmat(4:6,1:3);
P22 = Pmat(4:6,4:6);


ampr(1:3,ifre) = -P11\(P12*A_dw);
ampr(4:6,ifre) = A_dw;
ampt(1:3,ifre) = 0;
ampt(4:6,ifre) = P21*ampr(1:3,ifre)+P22*A_dw;

ampr(:,nfft+2-ifre) = conj(ampr(:,ifre));
ampt(:,nfft+2-ifre) = conj(ampt(:,ifre));


end


%% Wavelet convolution
rk_f = 100; %Hz
[w,tw] = ricker(rk_f,dt);

plw_rt = fft(ampr)*df;
plw_tt = fft(ampt)*df;

for i = 1:6
amprt(i,:) = conv(plw_rt(i,:),w');
amptt(i,:) = conv(plw_tt(i,:),w');
end



t_vec = layer(nlayer).tv;



%% Plot the figure
tt = [0:nfft-1]*dt;

figure
subplot(3,1,1)
plot(tt,amprt(1,1:nfft));
legend('S Slow Reflection Amplitude')
subplot(3,1,2)
plot(tt,amprt(2,1:nfft));
legend('S Fast Reflection Amplitude')
subplot(3,1,3)
plot(tt,amprt(3,1:nfft));
legend('P Reflection Amplitude')

figure
subplot(3,1,1)
plot(tt,amptt(4,1:nfft));
legend('P Transmission Amplitude')
subplot(3,1,2)
plot(tt,amptt(5,1:nfft));
legend('S Fast Transmission Amplitude')
subplot(3,1,3)
plot(tt,amptt(6,1:nfft));
legend('S Slow Transmission Amplitude')


