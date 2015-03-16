function [w,P,TM] = neoHookean(F,lambda,mu)
%NEOHOOKEAN takes the deformation gradient as an input and returns the
%strain energy density w, first Piola-Kirchhoff stress tensor P and
%tangent modulus TM as output.User must also specify the values of
%material constants 'lambda' and 'mu'

J = det(F);
Fi = inv(F);
C = F'*F;
I1 = C(1,1) + C(2,2) + C(3,3);

%Calculating strain energy density
w = (lambda/2)*(log(J))^2 - mu*(log(J)) + (mu/2)*(I1-3);

%Calculating the first Piola-Kirchoff Stress tensor
P = (lambda*log(J)- mu)*Fi' + mu*F;
% if(~isreal(P))
%     error('P has become imaginary.');
% end

%Assming 3-dimensions for a fourth-order tangent modulus
%TM = zeros(3,3,3,3);

%Calculating the TM components
% for l=1:3
%     for M=1:3
%         for n=1:3
%             for O=1:3
%                 TM(l,M,n,O) = lambda*Fi(M,l)*Fi(O,n) + mu*...
%                     (l==n)*(M==O) - (lambda*log(J)...
%                     - mu)*Fi(M,n)*Fi(O,l);
%             end
%         end
%     end
% end

%TM(l,M,n,O) = lambda*Fi(M,l)*Fi(O,n) + mu*...
%                     (l==n)*(M==O) - (lambda*log(J)...
%                     - mu)*Fi(M,n)*Fi(O,l);

V1 = permute(Fi,[2,1,3,4]);
V2 = permute(Fi,[3,4,2,1]);
V3 = permute(eye(3),[1,3,2,4]);
V4 = permute(eye(3),[3,1,4,2]);
V5 = permute(Fi,[3,1,2,4]);
V6 = permute(Fi,[2,3,4,1]);

TM = lambda*bsxfun(@times,V1,V2) + mu*bsxfun(@times,V3,V4) -...
    (lambda*log(J)- mu)*bsxfun(@times,V5,V6);


end