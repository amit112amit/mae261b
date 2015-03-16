function [w,S,C_IJKL] = neoHookeanSC(C,lambda,mu)
%NEOHOOKEAN takes the deformation gradient as an input and returns the
%strain energy density w, first Piola-Kirchhoff stress tensor P and
%tangent modulus TM as output.User must also specify the values of 
%material constants 'lambda' and 'mu'

I1 = C(1,1) + C(2,2) + C(3,3);
detF = sqrt(det(C));

dim = 3;

Identity = eye(3);
Ci = C\Identity; % Just taking inverse of C
%Cit = Ci.'; %C is symmetric so transpose(Ci) is same as Ci


%Calculating strain energy density
w = (lambda/2)*(log(detF))^2 - mu*(log(detF)) + (mu/2)*(I1-3);

%Calculating the second Piola-Kirchoff Stress tensor
S = (lambda*log(detF) - mu)*Ci + mu*Identity;

% Calculatiing C_IJKL by alternative analytical expression in terms of 
% w(C)
% We need to store the derivative of inverse transpose of C with C in Z
Z = zeros(dim,dim,dim,dim);
for i=1:3
    for j=1:3
        ei = Identity(:,i);
        ej = Identity(:,j);        
        Z(:,:,i,j) = (-Ci*ei)*((C\ej).'); % From matrix cook book
    end
end

C_IJKL = zeros(3,3,3,3);
for I=1:3
    for J=1:3
        for K=1:3
            for L=1:3  
               C_IJKL(I,J,K,L) = (lambda*log(detF)-mu)*Z(I,J,K,L) +...
                   lambda*Ci(I,J)*Ci(K,L)/2;
            end
        end
    end
end

end