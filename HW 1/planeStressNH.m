function [w,P,C,F33] = planeStressNH(D,lambda,mu,guessF33)
%PLANESTRESS takes F_{\alpha\beta} as input and returns strain energy w,
%  first Piola-Kirchhoff stress tensor P and tangent modulus C as output
%  consistent with plane stress condition

%Initial guess for F33 i.e. gamma and initialize P
gamma = guessF33;
tol = 10^-8;
iterLimit =100;
iterCount = 1;

%Newton iteration to solve for F33 such that P33 = 0
while(1)
    F = [D(1,1),D(1,2),0; D(2,1),D(2,2),0; 0,0,gamma];
    if(det(F) < 0)
        error('Invalid F. det(F) = %17.16f Iterations = %d.',...
            det(F),iterCount);
    end
    [~,P,TM] = neoHookean(F,lambda,mu);
    dgamma = -(P(3,3))/(TM(3,3,3,3));
    gamma = gamma + dgamma;
    iterCount = iterCount + 1;
%     if(~isreal(P))
%         error('P has become imaginary. Calculation terminated.');
%     end
    if(abs(P(3,3)) < tol || iterCount>iterLimit || abs(dgamma) < eps*10^3)
%         fprintf(['Newton iteration terminated.\n'...
%             'abs(P33) = %4.3e abs(dgamma) = %4.3e iterCount = %d\n'],...
%             abs(P(3,3)),abs(dgamma),iterCount);
        break;
    end
end
F33 = F(3,3);
[w,P,TM] = neoHookean(F,lambda,mu);
P = P(1:2,1:2);

%Corrections to the tangent modulus
C = TM(1:2,1:2,1:2,1:2);
for a=1:2
    for b=1:2
        for c=1:2
            for d=1:2
                C(a,b,c,d) = C(a,b,c,d) -...
                    TM(a,b,3,3)*TM(3,3,c,d)/TM(3,3,3,3);
            end
        end
    end
end

end