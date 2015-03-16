function [W,f,K] = T6QuadEle(X,x,quadOrder,lambda,mu)
% T6QUADELE takes as input X_{ia} and x_{ia}, the reference and spatial
% positions for each node of a triangular 6-node quadratic element and
% calculates strain energy W, internal nodal force array f and stiffness
% K. Input 'quadOrder' can be '1' or '2' for linear and quadratic Gauss
% quadrature respectively. lambda and mu are material properties.

dim = 2;
numNodes = 6;
guessF33 = 0.5;

% Get the quadrature points and weights
[points,weights] = TriGaussQuad(quadOrder);

W = 0;
f = zeros(dim,numNodes);
K = zeros(dim,numNodes,dim,numNodes);

for index = 1:length(weights)
    [~,DN] = T6quad(points(index,:));
    Jacobian = X*DN;
    DNDX = DN/Jacobian;
    
    % Calculate F_{iJ}(\theta_{index})
    F = x*DNDX;
    
    if(det(F) < 0)
        error('Deformation gradient has negative determinant.');
    end
    
    % Calculate strain energy
    [wq,P,~,~] = planeStressNH(F,lambda,mu,guessF33);
    W = W + wq*weights(index)*det(Jacobian);
    
    f = f + P*DNDX'*det(Jacobian)*weights(index);
end

for i=1:dim
    for a=1:numNodes
        for k=1:dim
            for b=1:numNodes
                
                for p=1:length(weights)
                    [~,DN] = T6quad(points(p,:));
                    Jacobian = X*DN;
                    DNDX = DN/Jacobian;
                    F = x*DNDX;                    
                    [~,~,C,~] = planeStressNH(F,lambda,mu,guessF33);
                    for J=1:dim
                        for L=1:dim
                            K(i,a,k,b) = K(i,a,k,b) +...
                                C(i,J,k,L)*DNDX(a,J)*...
                                DNDX(b,L)*det(Jacobian)*...
                                weights(p);
                        end
                    end
                end
                
            end
        end
    end
end


end