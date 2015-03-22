function [F,PKStress] = calcFandP(x,X,Lambda,numNodes,quadOrder,lame1,mu)
%calcFandP Calculates deformation gradient at quadrature point of a
%linear 3-node triangular element
%   x is the deformed position of nodes
%   X is the reference position of nodes
%   Lambda is the thickness stretch of the element
%   numNodes is the number of nodes in the element
%   quadOrder is order of quadrature
%
%   F is the deformation gradient in lab frame

dim = 3;

[points,weights] = TriGaussQuad(quadOrder);
F = zeros(dim,dim,numel(weights));
PKStress = zeros(dim,dim,numel(weights));

x = reshape(x,[dim,numNodes]);
X = reshape(X,[dim,numNodes]);

for z = 1:size(F,3)
    if numNodes == 3
        [~,DN] = T3Lin(points(:,z));
    elseif numNodes == 6
        [~,DN] = T6quad(points(:,z));
    end
    
    a_alpha_sub = x*DN;
    A_alpha_sub = X*DN;
    
    a_3 = cross(a_alpha_sub(:,1),a_alpha_sub(:,2));
    a_3 = a_3/norm(a_3);
    
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
    
    A_alpha_sup = G_dual(:,1:2);
    
    F(:,:,z) = a_alpha_sub(:,1)*(A_alpha_sup(:,1)).' +...
        a_alpha_sub(:,2)*(A_alpha_sup(:,2)).' + Lambda*a_3*A3.';
    [~,PKStress(:,:,z),~] = neoHookean(F(:,:,z),lame1,mu);
end

end

