% Compare for-loop result and bsxfun result
clear;
clc;

F = rand(3);

C_for  = zeros(3,3,3,3);
tic;
for l=1:3
    for M=1:3
        for n=1:3
            for O=1:3
                C_for(l,M,n,O) = F(M,l)*F(O,n);
            end
        end
    end
end
toc;
tic;

A = permute(F,[2,1,3,4]);
B = permute(F,[3,4,2,1]);
C_vec = bsxfun(@times,A,B);

toc;
for l=1:3
    for M=1:3
        for n=1:3
            for O=1:3
                fprintf('C_for(l,M,n,O):%d C_vec(l,M,n,O): %d\n',...
                    C_for(l,M,n,O),C_vec(l,M,n,O));
                if(C_for(l,M,n,O) ~= C_vec(l,M,n,O))
                    error('Alas!');
                end
            end
        end
    end
end