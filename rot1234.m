function cr = rot1234 (c, R)
% c is rank-4 stiffness tensor 
% cr is the tensor after rotation 
% yingcai zheng, 2/12/2019 
cr = c * 0;
for I=1:3
    for J=1: 3
        for K=1:3
            for L=1:3
%                 cr(I,J,K,L) =0;
                for i=1:3
                    for j=1:3
                        for k=1:3
                            for l=1:3
                                cr(I,J,K,L) = cr(I,J,K,L) ...
                                    + c(i,j,k,l)*R(I,i)*R(J,j)*R(K,k)*R(L,l);
                                
                            end
                        end
                    end
                end
                
            end
        end
    end
end
end
