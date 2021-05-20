%==============================================================================
% Compute Cij's from Thomsen's VTI parameters IN CWPmatlab
%==============================================================================

function [Cij] = thomsen2Cij_CWP(VP0, VS0, epsilon, delta, gamma,rho)

c33 = VP0^2*rho;   
c44 = VS0^2*rho;
c11 = c33*(1 + 2*epsilon);   
c13 = sqrt( 2*c33*(c33-c44)*delta + (c33-c44)^2 ) - c44;
c66 = c44*(1 + 2*gamma);

%==============================================================================
% Construct Cij's 
Cij = zeros(6);
Cij(1,1) = c11;   Cij(1,2) = c11 - 2*c66;   Cij(1,3) = c13;
                  Cij(2,2) = c11;           Cij(2,3) = c13;
                                            Cij(3,3) = c33;
Cij(4,4) = c44;   Cij(5,5) = c44;           Cij(6,6) = c66;

% Fill the symmetric part
for i=1:6
  for j=1:i-1
    Cij(i,j) = Cij(j,i);
  end;
end;

%==============================================================================

