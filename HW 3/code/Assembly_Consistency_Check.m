%% Check Assembly for Consistency

close all; clear; clc;
tic;

%************************** Generate the mesh ****************************%

% For a rectangular mesh
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

% Calculate global DOF number using global node number and local dof
% number.

a = BC_node_localDOF(:,1);
i = BC_node_localDOF(:,2); % The local DOF number at the global node


dofPerNode = 3;
% Convert the BC array in terms of the global DOF number for convenience.
prescribedDOF = [dofPerNode*(a-1)+i, BC_node_localDOF(:,3)];


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
f = [0;0;1];

% Quadrature order
quadOrder = 1;

% Thickness of the membrane
H = 1;

% Elastic constants
lambda = 5*10^8;
mu = 1.5*10^6;

% Set the initial guess for displacement
rng(0);
u = 0.001*rand(size(X,1),size(X,2));
%u = zeros(size(X,1),size(X,2)) % Uncomment for zero deformation case

% Using linear indexing as the global node number is same as linear index
u(prescribedDOF(:,1)) = prescribedDOF(:,2);

X = reshape(X,[dofPerNode*numNodes,1]);
u = reshape(u,[dofPerNode*numNodes,1]);

x = X + u;

% The exact values
[W_exact,r_exact,K_exact] = assemblyT3Lin(X,x,H,f,quadOrder,lambda,mu,...
    IEN,ID);

h = norm(x)*logspace(-7,-2);
errR = zeros(1,length(h));

% We will use the fact that
for z=1:length(h)
    r_approx = zeros(numel(unknownDOFs),1);
    for i=1:numel(unknownDOFs)
        x(unknownDOFs(i)) = x(unknownDOFs(i)) + h(z);
        [Wp,~,~] = assemblyT3Lin(X,x,H,f,quadOrder,lambda,mu,IEN,ID);
        x(unknownDOFs(i)) = x(unknownDOFs(i)) - 2*h(z);
        [Wm,~,~] = assemblyT3Lin(X,x,H,f,quadOrder,lambda,mu,IEN,ID);
        x(unknownDOFs(i)) = x(unknownDOFs(i)) + h(z);
        r_approx(i) = (Wp - Wm)/(2*h(z));
    end
    errR(z) = max(max(abs(r_exact-r_approx)))/norm(r_exact);
    %errR(z) = max(max(abs(r_exact-r_approx)));
end

slope1 = (log(errR(50))-log(errR(45)))/(log(h(50))-log(h(45)));
fprintf('Slope for error in f_internal: %17.16f\n',slope1);

errK = zeros(1,length(h));
for z=1:length(h)
    K_approx = zeros(numel(unknownDOFs),numel(unknownDOFs));
    for m=1:numel(unknownDOFs)
        for n=1:numel(unknownDOFs)      
            x(unknownDOFs(n)) = x(unknownDOFs(n)) + h(z);
            [~,rp,~] = assemblyT3Lin(X,x,H,f,quadOrder,...
                lambda,mu,IEN,ID);
            x(unknownDOFs(n)) = x(unknownDOFs(n)) - 2*h(z);
            [~,rm,~] = assemblyT3Lin(X,x,H,f,quadOrder,...
                lambda,mu,IEN,ID);
            x(unknownDOFs(n)) = x(unknownDOFs(n)) + h(z);
            K_approx(m,n) = (rp(m) - rm(m))/(2*h(z));
        end
    end
    errK(z) = max(max(max(max(abs((K_exact - K_approx))))))/...
        max(max(max(max(K_exact))));
    %errK(z) = max(max(max(max(abs((K_exact - K_approx))))));
end
figure(2);
loglog(h,errR,h,errK);
xlabel('log(h)');
ylabel('log(error)')
legend('Log of error in f_{int}','Log of error in K_{iakb}');
slope2 = (log(errK(50))-log(errK(45)))/(log(h(50))-log(h(45)));
fprintf('Slope for error in K: %17.16f\n',slope2);

toc;
