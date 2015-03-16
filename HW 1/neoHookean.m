function [w,P,TM] = neoHookean(F,lambda,mu)
%NEOHOOKEAN takes the deformation gradient as an input and returns the
%strain energy density w, first Piola-Kirchhoff stress tensor P and
%tangent modulus TM as output.User must also specify the values of 
%material constants 'lambda' and 'mu'

J = det(F);
Fi = inv(F);
C = F'*F;
I1 = C(1,1) + C(2,2) + C(3,3);

%kronDel = @(x,y)(x==y);

%Calculating strain energy density
w = (lambda/2)*(log(J))^2 - mu*(log(J)) + (mu/2)*(I1-3);

%Calculating the first Piola-Kirchoff Stress tensor
P = (lambda*log(J)- mu)*Fi' + mu*F;

%Assming 3-dimensions for a fourth-order tangent modulus
TM = zeros(3,3,3,3);

%Calculating the TM components
for l=1:3
    for M=1:3
        for n=1:3
            for O=1:3
                TM(l,M,n,O) = lambda*Fi(M,l)*Fi(O,n) + mu*...
                    (l==n)*(M==O) - (lambda*log(J)...
                    - mu)*Fi(M,n)*Fi(O,l);
            end
        end
    end
end
end