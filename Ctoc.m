function c = Ctoc(C)
%% transfer C_ij to c_ijkl 
%% input C_ij
%% output c_ijkl

for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                % I
                if i==1 && j==1, I = 1; end
                if i==2 && j==2, I = 2; end
                if i==3 && j==3, I = 3; end
                if i==2 && j==3, I = 4; end
                if i==3 && j==2, I = 4; end  % sym
                if i==1 && j==3, I = 5; end
                if i==3 && j==1, I = 5; end  % sym
                if i==1 && j==2, I = 6; end
                if i==2 && j==1, I = 6; end  % sym
                % J
                if k==1 && l==1, J = 1; end
                if k==2 && l==2, J = 2; end
                if k==3 && l==3, J = 3; end
                if k==2 && l==3, J = 4; end
                if k==3 && l==2, J = 4; end  % sym
                if k==1 && l==3, J = 5; end
                if k==3 && l==1, J = 5; end  % sym
                if k==1 && l==2, J = 6; end
                if k==2 && l==1, J = 6; end  % sym
                c(i,j,k,l) = C(I, J);
            end
        end
    end
end




end