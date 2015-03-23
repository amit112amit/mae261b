function [StrainEnergy,fint_ia,fext_ia,Kiakb,L] = T3MembraneEle(X,x,H,f,...
    quadOrder,lambda,mu,thicknessStretch)
%T3MEMBRANEELE Calculates strain energy, internal force vector and
%stiffness modulus for the triangular 3 node membrane element
%   X is reference configuration, x is deformed configuration, H is
%   thickness, f is constant distributed load
%   W, fia, Kiakb are strain energy, internal force and stiffness modulus
%   respectively
%   L is the thickness stretch, Lambda. (not the Lame coeff 'l'ambda)

dim = 3;
numNodes = 3;

x = reshape(x,[dim,numNodes]);
X = reshape(X,[dim,numNodes]);

isPlaneStress = false;

% Get the quadrature points and weights
[points,weights] = TriGaussQuad(quadOrder);

StrainEnergy = 0;
fint_ia = zeros(dim,numNodes);
fext_ia = zeros(dim,numNodes);
Kiakb = zeros(dim,numNodes,dim,numNodes);

for z=1:length(weights)
    % Get the shape functions and their derivatives for T3Lin
    [N,DN] = T3Lin(points(:,z));
    a_alpha_sub = x*DN;
    A_alpha_sub = X*DN;
    
    % Calculate thickness stretch L
    [strEngDen,~,Ctilda,otherData]= calcAllCs(a_alpha_sub,...
        A_alpha_sub,thicknessStretch,H,lambda,mu,isPlaneStress);
    sqrt_A = otherData.sqrt_A;
    
    % Calculate the strain energy from the density
    StrainEnergy = StrainEnergy +...
        strEngDen*sqrt_A*weights(z)*H;
    
    % Calculate internal forces
    n_alpha = otherData.n_alpha;
    fint_ia = fint_ia + n_alpha*(DN.')*sqrt_A*weights(z);
    
    % Calculate external forces
    fext_ia = fext_ia + f*(N.')*sqrt_A*weights(z);
    
    %Calculate K_material
    tau = otherData.tau;
    %fprintf('Tau(3,3) = %17.16f\n',tau(3,3));
    
    K_material = zeros(dim,numNodes,dim,numNodes);
    
    for i=1:dim
        for a=1:numNodes
            for k=1:dim
                for b=1:numNodes
                    
                    for p=1:2
                        for q=1:2
                            for r=1:2
                                for s=1:2
                                    
                                    tmp = a_alpha_sub(:,q)*...
                                        (a_alpha_sub(:,s).');
                                    
                                    K_material(i,a,k,b) =...
                                        K_material(i,a,k,b) + ...
                                        2*Ctilda(p,q,r,s)*tmp(i,k)*...
                                        DN(a,p)*DN(b,r);
                                end
                            end
                        end
                    end
                    
                end
            end
        end
    end
    
    K_geometric = zeros(dim,numNodes,dim,numNodes);
    
    for i=1:dim
        for a=1:numNodes
            for k=1:dim
                for b=1:numNodes
                    
                    for p=1:2
                        for q=1:2
                            K_geometric(i,a,k,b) = K_geometric(i,a,k,b)...
                                + tau(p,q)*(i==k)*DN(a,p)*DN(b,q);                            
                        end
                    end
                    
                end
            end
        end
    end
    
    Kiakb = Kiakb + (K_material + K_geometric)*sqrt_A*H*weights(z);
    
end

% Reshape all the output matrices
Kiakb = reshape(Kiakb,[dim*numNodes,dim*numNodes]);
fint_ia = reshape(fint_ia,[dim*numNodes,1]);
fext_ia = reshape(fext_ia,[dim*numNodes,1]);
L = otherData.Lambda;
end

