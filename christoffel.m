function [p_vec, t_vec] = christoffel(rho, Cij, px, py)
%% Solve the Green-Christoffel equation
%% Input: rho density of the layer
%% Input: Cij the 2-order elastic tensor
%% Input: px, py, in-plane wavenumber (given when initialize the plane wave)
%% Output: p_vec, 6*3 matrix, 6 for 3 upgoing 3 downgoing waves.
%%         3 for px py pz components
%% Output: t_vec, polarization matrix, each row is a polarization vector
%%         6*3 matrix, meaning of rwo and column is simular with p_vec.

%% First, calculate cijkl*pj*pl
cijkl = C2toc4(Cij); % transfer Cij to cijkl

[Q,P,M] = get_coeffmat(cijkl, px, py); % Q coef p^0 P coef p^1 M coef p^2

m11 = M(1,1); m12 = M(1,2); m13=M(1,3);
m21 = M(2,1); m22 = M(2,2); m23=M(2,3);
m31 = M(3,1); m32 = M(3,2); m33=M(3,3);

p11 = P(1,1); p12 = P(1,2); p13=P(1,3);
p21 = P(2,1); p22 = P(2,2); p23=P(2,3);
p31 = P(3,1); p32 = P(3,2); p33=P(3,3);

q11 = Q(1,1); q12 = Q(1,2); q13=Q(1,3);
q21 = Q(2,1); q22 = Q(2,2); q23=Q(2,3);
q31 = Q(3,1); q32 = Q(3,2); q33=Q(3,3);

pw = rho;

pw_m = [pw,0,0;0,pw,0;0,0,pw];

c6 = (-1).*m13.*m22.*m31+m12.*m23.*m31+m13.*m21.*m32+(-1).*m11.*m23.* ...
    m32+(-1).*m12.*m21.*m33+m11.*m22.*m33;

c5 = (-1).*m23.*m32.*p11+m22.*m33.*p11+m23.*m31.*p12+(-1).*m21.*m33.* ...
    p12+(-1).*m22.*m31.*p13+m21.*m32.*p13+m13.*m32.*p21+(-1).*m12.* ...
    m33.*p21+(-1).*m13.*m31.*p22+m11.*m33.*p22+m12.*m31.*p23+(-1).* ...
    m11.*m32.*p23+(-1).*m13.*m22.*p31+m12.*m23.*p31+m13.*m21.*p32+(-1) ...
    .*m11.*m23.*p32+(-1).*m12.*m21.*p33+m11.*m22.*p33;

c4 = (-1).*m33.*p12.*p21+m32.*p13.*p21+m33.*p11.*p22+(-1).*m31.*p13.* ...
    p22+(-1).*m32.*p11.*p23+m31.*p12.*p23+m23.*p12.*p31+(-1).*m22.* ...
    p13.*p31+(-1).*m13.*p22.*p31+m12.*p23.*p31+(-1).*m23.*p11.*p32+ ...
    m21.*p13.*p32+m13.*p21.*p32+(-1).*m11.*p23.*p32+m22.*p11.*p33+(-1) ...
    .*m21.*p12.*p33+(-1).*m12.*p21.*p33+m11.*p22.*p33+m12.*m21.*pw+( ...
    -1).*m11.*m22.*pw+m13.*m31.*pw+m23.*m32.*pw+(-1).*m11.*m33.*pw+( ...
    -1).*m22.*m33.*pw+(-1).*m23.*m32.*q11+m22.*m33.*q11+m23.*m31.*q12+ ...
    (-1).*m21.*m33.*q12+(-1).*m22.*m31.*q13+m21.*m32.*q13+m13.*m32.* ...
    q21+(-1).*m12.*m33.*q21+(-1).*m13.*m31.*q22+m11.*m33.*q22+m12.* ...
    m31.*q23+(-1).*m11.*m32.*q23+(-1).*m13.*m22.*q31+m12.*m23.*q31+ ...
    m13.*m21.*q32+(-1).*m11.*m23.*q32+(-1).*m12.*m21.*q33+m11.*m22.* ...
    q33;

c3 = (-1).*p13.*p22.*p31+p12.*p23.*p31+p13.*p21.*p32+(-1).*p11.*p23.* ...
    p32+(-1).*p12.*p21.*p33+p11.*p22.*p33+(-1).*m22.*p11.*pw+(-1).* ...
    m33.*p11.*pw+m21.*p12.*pw+m31.*p13.*pw+m12.*p21.*pw+(-1).*m11.* ...
    p22.*pw+(-1).*m33.*p22.*pw+m32.*p23.*pw+m13.*p31.*pw+m23.*p32.*pw+ ...
    (-1).*m11.*p33.*pw+(-1).*m22.*p33.*pw+m33.*p22.*q11+(-1).*m32.* ...
    p23.*q11+(-1).*m23.*p32.*q11+m22.*p33.*q11+(-1).*m33.*p21.*q12+ ...
    m31.*p23.*q12+m23.*p31.*q12+(-1).*m21.*p33.*q12+m32.*p21.*q13+(-1) ...
    .*m31.*p22.*q13+(-1).*m22.*p31.*q13+m21.*p32.*q13+(-1).*m33.*p12.* ...
    q21+m32.*p13.*q21+m13.*p32.*q21+(-1).*m12.*p33.*q21+m33.*p11.*q22+ ...
    (-1).*m31.*p13.*q22+(-1).*m13.*p31.*q22+m11.*p33.*q22+(-1).*m32.* ...
    p11.*q23+m31.*p12.*q23+m12.*p31.*q23+(-1).*m11.*p32.*q23+m23.* ...
    p12.*q31+(-1).*m22.*p13.*q31+(-1).*m13.*p22.*q31+m12.*p23.*q31+( ...
    -1).*m23.*p11.*q32+m21.*p13.*q32+m13.*p21.*q32+(-1).*m11.*p23.* ...
    q32+m22.*p11.*q33+(-1).*m21.*p12.*q33+(-1).*m12.*p21.*q33+m11.* ...
    p22.*q33;

c2 = p12.*p21.*pw+(-1).*p11.*p22.*pw+p13.*p31.*pw+p23.*p32.*pw+(-1).* ...
    p11.*p33.*pw+(-1).*p22.*p33.*pw+m11.*pw.^2+m22.*pw.^2+m33.*pw.^2+( ...
    -1).*p23.*p32.*q11+p22.*p33.*q11+(-1).*m22.*pw.*q11+(-1).*m33.* ...
    pw.*q11+p23.*p31.*q12+(-1).*p21.*p33.*q12+m21.*pw.*q12+(-1).*p22.* ...
    p31.*q13+p21.*p32.*q13+m31.*pw.*q13+p13.*p32.*q21+(-1).*p12.*p33.* ...
    q21+m12.*pw.*q21+(-1).*m33.*q12.*q21+m32.*q13.*q21+(-1).*p13.* ...
    p31.*q22+p11.*p33.*q22+(-1).*m11.*pw.*q22+(-1).*m33.*pw.*q22+m33.* ...
    q11.*q22+(-1).*m31.*q13.*q22+p12.*p31.*q23+(-1).*p11.*p32.*q23+ ...
    m32.*pw.*q23+(-1).*m32.*q11.*q23+m31.*q12.*q23+(-1).*p13.*p22.* ...
    q31+p12.*p23.*q31+m13.*pw.*q31+m23.*q12.*q31+(-1).*m22.*q13.*q31+( ...
    -1).*m13.*q22.*q31+m12.*q23.*q31+p13.*p21.*q32+(-1).*p11.*p23.* ...
    q32+m23.*pw.*q32+(-1).*m23.*q11.*q32+m21.*q13.*q32+m13.*q21.*q32+( ...
    -1).*m11.*q23.*q32+(-1).*p12.*p21.*q33+p11.*p22.*q33+(-1).*m11.* ...
    pw.*q33+(-1).*m22.*pw.*q33+m22.*q11.*q33+(-1).*m21.*q12.*q33+(-1) ...
    .*m12.*q21.*q33+m11.*q22.*q33;

c1 = p11.*pw.^2+p22.*pw.^2+p33.*pw.^2+(-1).*p22.*pw.*q11+(-1).*p33.* ...
    pw.*q11+p21.*pw.*q12+p31.*pw.*q13+p12.*pw.*q21+(-1).*p33.*q12.* ...
    q21+p32.*q13.*q21+(-1).*p11.*pw.*q22+(-1).*p33.*pw.*q22+p33.*q11.* ...
    q22+(-1).*p31.*q13.*q22+p32.*pw.*q23+(-1).*p32.*q11.*q23+p31.* ...
    q12.*q23+p13.*pw.*q31+p23.*q12.*q31+(-1).*p22.*q13.*q31+(-1).* ...
    p13.*q22.*q31+p12.*q23.*q31+p23.*pw.*q32+(-1).*p23.*q11.*q32+p21.* ...
    q13.*q32+p13.*q21.*q32+(-1).*p11.*q23.*q32+(-1).*p11.*pw.*q33+(-1) ...
    .*p22.*pw.*q33+p22.*q11.*q33+(-1).*p21.*q12.*q33+(-1).*p12.*q21.* ...
    q33+p11.*q22.*q33;

c0 = (-1).*pw.^3+pw.^2.*q11+pw.*q12.*q21+pw.^2.*q22+(-1).*pw.*q11.*q22+ ...
    pw.*q13.*q31+(-1).*q13.*q22.*q31+q12.*q23.*q31+q13.*q21.*q32+pw.* ...
    q23.*q32+(-1).*q11.*q23.*q32+pw.^2.*q33+(-1).*pw.*q11.*q33+(-1).* ...
    q12.*q21.*q33+(-1).*pw.*q22.*q33+q11.*q22.*q33;

pol =[c6 c5 c4 c3 c2 c1 c0];

pz = transpose(roots (pol));

[~,indx_sort] = sort(real(pz));

pz = pz(indx_sort);



%% Debug interval 1
% p_t = pz(1);
% poly_p = [p_t^6,p_t^5,p_t^4,p_t^3,p_t^2,p_t,1];
% res = dot(poly_p,pol)
% GC = p_t^2*M + p_t*P + Q - pw_m;
% %det(GC)
% res - det(GC)

p_vec = zeros(6,3);
p_vec(:,1) = px;
p_vec(:,2) = py;
p_vec(:,3) = pz;
t_vec = GC_solve(M,P,Q,pw_m,pz);
    
    
    
end





