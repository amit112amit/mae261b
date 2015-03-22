function [W,r,kiakb,L] = assemblyT3Lin(X,x,H,f,quadOrder,lambda,mu,IEN,...
    ID,L)
%AssemblyT3Lin Performs assembly and returns potential energy, residual
%force and corresponding stiffness.
% IEN is the connectivity matrix of the mesh
% ID gives the global eqn number for every global DOF
% Other input parameters are same as for T3MembraneEle
% H should be an array for thickness of all elements ordered as per element
% number
% L is an array of thickness stretches for each element
%
% W is the potential energy
% r is the residual of the fint - fext
% kiakb is the stiffness of the whole system
% H is the updated thickness for all elements

numEle = size(IEN,1); % Number of rows in connectivity table
nodesPerEle = size(IEN,2);
dofPerNode = 3;

numUnknownDOF = max(ID(:,2));

kiakb = zeros(numUnknownDOF,numUnknownDOF);
fi_global = zeros(numUnknownDOF,1);
fe_global = zeros(numUnknownDOF,1);
W = 0;

tmp2 = (1:dofPerNode).'; % We will use it to calculate global DOF number
tmp2 = repmat(tmp2,[nodesPerEle,1]);

for i=1:numEle
    eleNodeNum = IEN(i,:);
    
    % Fancy stuff to calculate globalDOF numbers for the global node
    % numbers without using for-loop
    tmp1 = repmat(eleNodeNum,[dofPerNode,1]);
    tmp1 = reshape(tmp1,[numel(tmp1),1]);
    eleGlobalDOF = dofPerNode*(tmp1-1) + tmp2;
    %eleLocalDOF = (1:dofPerNode*nodesPerEle).';
    
    Xele = X(eleGlobalDOF);
    xele = x(eleGlobalDOF);
    uele = xele - Xele; % displacement    
    
    [W_ele_int,fi_ele,fe_ele,K_ele,L(i)] = T3MembraneEle(Xele,xele,H(i),...
        f,quadOrder,lambda,mu,L(i));    
    
    % Including the energy contribution from constrained DOFs too.
    W_ele_ext = fe_ele.'*uele;    
    W = W + W_ele_int - W_ele_ext;
    
    globalEqNums = ID(eleGlobalDOF,2);
    logicalIndex = globalEqNums > 0;
    globalEqNums = globalEqNums(logicalIndex);  
    
    kiakb(globalEqNums,globalEqNums) = kiakb(globalEqNums,...
        globalEqNums) + K_ele(logicalIndex,logicalIndex);
    
    fi_global(globalEqNums) = fi_global(globalEqNums) +...
        fi_ele(logicalIndex);
    
    fe_global(globalEqNums) = fe_global(globalEqNums) +...
        fe_ele(logicalIndex);  
    
end
r = fi_global - fe_global;

end

