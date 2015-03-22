%% Check rank of stiffness matrix after assembly

close all; clear; clc;
tic;

%************************** Generate the mesh ****************************%

% For rectangular mesh
xlim = 1;
ylim = 1;
zlim = 1;
[mesh_x,mesh_y,mesh_z] = meshgrid(0:xlim,0:ylim,1:zlim);

% For a single triangle
% mesh_x = [1;0;1];
% mesh_y = [1;0;0];
% mesh_z = [1;1;1];

% Generate the mesh
IEN = delaunay(mesh_x,mesh_y);

% Convert the mesh grid points into x_ia
numNodes = numel(mesh_x);

% The reference configuration
X = [reshape(mesh_x,numNodes,1), reshape(mesh_y,numNodes,1),...
    reshape(mesh_z,numNodes,1)].';

trisurf(IEN,X(1,:),X(2,:),X(3,:));


%****************** Boundary Condition and Meta Arrays *******************%

% Need to know the Global DOF number that has to be constrained.

% Boundary condition specified with help of global node number and the
% local dof number at that node. So local dof at that node goes from 1 to
% dofpernode. In the following array, first column is global node number,
% second column is local dof number, third column is the prescribed value.

BC_node_localDOF = [1, 1, 0; 1, 2, 0; 1, 3, 0];
%BC_node_localDOF = [1,1,0;1,2,0;2,1,0;2,2,0;3,1,0;3,2,0;];
%BC_node_localDOF = [];

% Calculate global DOF number using global node number and local dof
% number.

a = BC_node_localDOF(:,1);
i = BC_node_localDOF(:,2); % The local DOF number at the global node


dofPerNode = 3;
% Convert the BC array in terms of the global DOF number for convenience.
prescribedDOF = [dofPerNode*(a-1)+i, BC_node_localDOF(:,3)];
%prescribedDOF = [];

ID = zeros(numNodes*dofPerNode,2);
ID(:,1) = 1:numNodes*dofPerNode;

unknownDOFs = (1:numNodes*dofPerNode).';
unknownDOFs(prescribedDOF(:,1)) = [];

globalEqNum = 1;
for i = 1:numel(unknownDOFs)
    ID(unknownDOFs(i),2) = globalEqNum;
    globalEqNum = globalEqNum + 1;
end

%******************************** Assembly *******************************%

% External transverse load
f = [0;0;0];

% Quadrature order
quadOrder = 1;

% Thickness of the membrane
H = 1;

% Elastic constants
lambda = 5*10^8;
mu = 1.5*10^6;

% Set the initial guess for displacement
rng(0);

u = zeros(size(X,1),size(X,2)); % Uncomment for zero deformation case
u(3,:) = 0.001*rand(1,size(X,2));% To avoid zero stiffness in z-direction

% Using linear indexing as the global node number is same as linear index
%u(prescribedDOF(:,1)) = prescribedDOF(:,2);

X = reshape(X,[dofPerNode*numNodes,1]);
u = reshape(u,[dofPerNode*numNodes,1]);

x = X + u;

%
[~,~,K1] = assemblyT3Lin(X,x,H,f,quadOrder,lambda,mu,...
    IEN,ID);
r1 = findStiffnessRank(K1);
fprintf('Zero deformation Stiffness Matrix rank = %d\n\n',r1);

u = 0.001*rand(size(X,1),size(X,2));
x = X + u;

[~,~,K2] = assemblyT3Lin(X,x,H,f,quadOrder,lambda,mu,...
    IEN,ID);

r2 = findStiffnessRank(K2);
fprintf('Finite deformation Stiffness Matrix rank = %d\n',r2);