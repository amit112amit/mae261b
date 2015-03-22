%% Explore delaunay()
close all; clear; clc;
xlim = 1;
ylim = 1;
zlim = 1;
[mesh_x,mesh_y,mesh_z] = meshgrid(0:xlim,0:ylim,1:zlim);
IEN = delaunay(mesh_x,mesh_y);

% Convert the mesh grid points into x_ia
numNodes = numel(mesh_x);
x_ia = [reshape(mesh_x,numNodes,1), reshape(mesh_y,numNodes,1),...
    reshape(mesh_z,numNodes,1)].';

% Checking if the Global DOF number is same as the linear index.
for i=1:3
    for j=1:numNodes
        tmp = sub2ind(size(x_ia),i,j);
        fprintf('Global DOF: %d , Linear Index: %d\n',3*(j-1)+i,tmp);
    end
end

trisurf(IEN,x_ia(1,:),x_ia(2,:),x_ia(3,:));
% Conclusion: TRI is exactly the matrix IEN that Hughes wants

%% How to construct ID matrix ?
% This matrix stores the global equation number for every global degree of
% freedom
% ID(a) =  P if this is not a prescribed DOF otherwise its 0 where i is the
% global degree of freedom number.
%
% Global node number a will have global DOF numbers 3*(a-1)+i where i=1,2,3
%
% Input: We need to know which global DOFs are prescribed.


% We will store the prescribed DOFs in an array of numBC-by-2. First column
% indicates the global dof number, second column contains the prescribed
% value of the dof.
numBC = 2;
prescribedDOF = [1,0;...
                2,0;...
                3,0];

dofPerNode = 3;

ID = zeros(numNodes*dofPerNode,2);
ID(:,1) = 1:numNodes*dofPerNode;

unknownDOFs = (1:numNodes*dofPerNode).';
unknownDOFs(prescribedDOF(:,1)) = [];

globalEqNum = 1;
for i = 1:numel(unknownDOFs)
    ID(unknownDOFs(i),2) = globalEqNum;
    globalEqNum = globalEqNum + 1;
end


