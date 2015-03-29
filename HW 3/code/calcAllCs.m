function [strEngDen,PKstress,Ctilda,otherData] = calcAllCs(a_alpha_sub,...
    A_alpha_sub,thicknessStretch,H,lame1,mu,isPlaneStress)
%CALCALLCS Planes stress neo-Hookean for curvilinear co-ordinates
%   a_alpha is tangent basis vector for mid-plane in spatial configuration
%   A_alpha is dual basis vector for mid-plane in reference configuration
%   thicknessStretch is thickness stretch
%   strEngDen is strain energy density and
%   PKstress is Piola-Kirchhoff Stress

dim = 3;

if(~isPlaneStress)
    Lambda = 1;
else
    Lambda = thicknessStretch; % The lambda from our notation in notes
end

% a_3 is same as a3
a_3 = cross(a_alpha_sub(:,1),a_alpha_sub(:,2));
a_3 = a_3/norm(a_3);

g = [a_alpha_sub,a_3];

g_ij = g.'*g; % metric tensor

gij = inv(g_ij); % dual metric tensor

g_dual = zeros(3);
for i=1:3
    g_dual(:,i) = gij(i,1)*g(:,1) + gij(i,2)*g(:,2) + gij(i,3)*g(:,3);
end

% A3 is same as A_3
A3 = cross(A_alpha_sub(:,1),A_alpha_sub(:,2));
A3 = A3/norm(A3);

G = [A_alpha_sub,A3];

G_ij = G.'*G; % metric tensor

Gij = inv(G_ij); % dual metric tensor
G_dual = zeros(3);
for i=1:3
    G_dual(:,i) = Gij(i,1)*G(:,1) + Gij(i,2)*G(:,2) + Gij(i,3)*G(:,3);
end

sqrt_A = sqrt(det(G_ij));

A_alpha_sup = G_dual(:,1:2);


if(isPlaneStress)
    maxIter = 100;
    iterCount = 0;
    tol = mu*10^(-15);
    
    retry = 0;
    while(iterCount < maxIter)
        % Calculate deformation gradient
        F = a_alpha_sub(:,1)*(A_alpha_sup(:,1)).' +...
            a_alpha_sub(:,2)*(A_alpha_sup(:,2)).' + Lambda*a_3*A3.';
        
        [strEngDen,PKstress,C_iJkL] = neoHookean(F,lame1,mu);
        
        T = a_3'*(PKstress*A3);
        
        if(abs(T) < tol)
            %fprintf('Tolerance met!\n');
            break;
        end
        
        % Calculate component of C_iJkL in the curvilinear frame
        C_3333 = 0;
        for I=1:dim
            for J=1:dim
                for K=1:dim
                    for L=1:dim
                        C_3333 = C_3333 + a_3(I)*A3(J)*C_iJkL(I,J,K,L)*...
                            a_3(K)*A3(L);
                    end
                end
            end
        end
        
        
        % The Newton iteration update
        dLambda = -T/C_3333;
        
        if(abs(dLambda) < eps)
            %fprintf('Newton update is too small!\n');
            break;
        end
        Lambda = Lambda + dLambda;
        
        % Don't let Lambda become negative!!!
        if (Lambda < 0)
            fprintf('calcAllCs(): Retry %d\n',retry+1);
            switch retry
                case 0
                    Lambda = 0.005;
                    retry = 1;
                case 1
                    Lambda = 0.00005;
                    retry = 2;
                case 2
                    Lambda = 0.00000005;
                    retry = 3;
                otherwise
                    fprintf('Failed after 3 attemtpts! Please debug.\n');
            end
        end
        iterCount = iterCount + 1;
    end
    
    if(iterCount >= maxIter && abs(T)>100000*tol)
        fprintf(['**** Plane stress Newton ',...
            'iterations did not converge. ****\n T = %17.16f',...
            ' dLambda = %17.16f\n'],T,dLambda);
    end
else
    F = a_alpha_sub(:,1)*(A_alpha_sup(:,1)).' +...
        a_alpha_sub(:,2)*(A_alpha_sup(:,2)).' + Lambda*a_3*A3.';
    
    [strEngDen,PKstress,C_iJkL] = neoHookean(F,lame1,mu);
    
end

% To derive C_IJKL from C_iJkL
C_IJKL = zeros(dim,dim,dim,dim);
Finv = inv(F);
S = F\PKstress;
for I=1:dim
    for J=1:dim
        for K=1:dim
            for L=1:dim
                for p=1:dim
                    for q=1:dim
                        C_IJKL(I,J,K,L) = C_IJKL(I,J,K,L) +...
                            0.5*Finv(I,p)*Finv(K,q)*(C_iJkL(p,J,q,L)-...
                            (p==q)*S(J,L));
                    end
                end
            end
        end
    end
end

Cijkl = zeros(dim,dim,dim,dim);
for i=1:dim
    for j=1:dim
        for k=1:dim
            for l=1:dim
                
                for I=1:dim
                    for J=1:dim
                        for K=1:dim
                            for L=1:dim
                                Cijkl(i,j,k,l) = Cijkl(i,j,k,l) +...
                                    C_IJKL(I,J,K,L)*...
                                    G_dual(I,i)*G_dual(J,j)*...
                                    G_dual(K,k)*G_dual(L,l);
                            end
                        end
                    end
                end
                
            end
        end
    end
end

if(isPlaneStress)
    Ctilda = zeros(2,2,2,2);
    for p=1:2
        for q=1:2
            for r=1:2
                for s=1:2
                    Ctilda(p,q,r,s) = Cijkl(p,q,r,s) - Cijkl(p,q,3,3)*...
                        Cijkl(3,3,r,s)/Cijkl(3,3,3,3);
                end
            end
        end
    end
else
    Ctilda = Cijkl(1:2,1:2,1:2,1:2);
end

% Calculating the components of Tau in curvilinear frame
tauij = zeros(3);
for i=1:3
    for j=1:3
        tmp = S*G_dual(:,j);
        tauij(i,j) = dot(G_dual(:,i),tmp);
    end
end

otherData.sqrt_A = sqrt_A;
otherData.tau = tauij;
otherData.n_alpha = PKstress*A_alpha_sup*H;
otherData.Lambda = Lambda;

end

