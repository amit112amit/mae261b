%% Equibi-axial strain

close all; clear; clc;
tic;

%************************** Generate the mesh ****************************%

% For rectangular mesh
xlim = 3;
ylim = 3;
zlim = 1;
[mesh_x,mesh_y,mesh_z] = meshgrid(0:xlim,0:ylim,1:zlim);

% Generate the mesh
IEN = delaunay(mesh_x,mesh_y);

% Convert the mesh grid points into x_ia
numNodes = numel(mesh_x);

% The reference configuration
X = [reshape(mesh_x,numNodes,1), reshape(mesh_y,numNodes,1),...
    reshape(mesh_z,numNodes,1)].';

%trisurf(IEN,X(1,:),X(2,:),X(3,:));


%****************** Boundary Condition and Meta Arrays *******************%

% Need to know the Global DOF number that has to be constrained. We will
% also constrain the plate to remain stationary in the z-direction

dofPerNode = 3;

constrainedSide1 = zeros(3,xlim+1);
constrainedSide1(1,:) = 0:xlim;
constrainedSide1(3,:) = ones(1,xlim+1);

constrainedSide2 = zeros(3,ylim+1);
constrainedSide2(1,:) = xlim*ones(1,xlim+1);
constrainedSide2(2,:) = 0:ylim;
constrainedSide2(3,:) = ones(1,xlim+1);

constrainedSide3 = zeros(3,xlim+1);
constrainedSide3(1,:) = 0:xlim;
constrainedSide3(2,:) = ylim*ones(1,ylim+1);
constrainedSide3(3,:) = ones(1,xlim+1);

constrainedSide4 = zeros(3,ylim+1);
constrainedSide4(1,:) = zeros(1,ylim+1);
constrainedSide4(2,:) = 0:ylim;
constrainedSide4(3,:) = ones(1,xlim+1);

BC = cell(1,3);
BC{1,1} = [constrainedSide1.';constrainedSide2.';constrainedSide3.';...
    constrainedSide4.';X.'];

BC{1,2} = [2*ones(xlim+1,1);ones(ylim+1,1);2*ones(xlim+1,1);...
    ones(ylim+1,1);3*ones(size(X,2),1)];

u_BC = (0.5)*xlim;

BC{1,3} = [-1*u_BC*ones(xlim+1,1);u_BC*ones(ylim+1,1);...
    u_BC*ones(xlim+1,1);-1*u_BC*ones(ylim+1,1);zeros(size(X,2),1)];

prescribedDOF = getBCmatrix(BC,X);

% Generate ID matrix
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
H = 0.001*xlim*ones(size(IEN,1),1);
thicknessStretch = ones(size(IEN,1),1);

% Elastic constants
lambda = 5*10^8;
mu = 1.5*10^6;

% Set the initial guess for displacement
rng(0);

u = zeros(size(X,1),size(X,2)); % Uncomment for zero deformation case

% Using linear indexing as the global node number is same as linear index
u(prescribedDOF(:,1)) = prescribedDOF(:,2);

X = reshape(X,[dofPerNode*numNodes,1]);
u_total = reshape(u,[dofPerNode*numNodes,1]);

u_steps = 10;
u = u_total/u_steps;

tol = norm(mu)*10^(-15);

X_orig = X; % Back up original reference config
H_orig = H;

for i=1:u_steps
    x = X + u;   
    maxIter = 100;    
    % Newton Iterations   
    while(1)
        [W,r,kiakb,L] = assemblyT3Lin(X,x,H,f,quadOrder,lambda,mu,...
            IEN,ID);
        kiakb = sparse(kiakb);
        u_Newton = -kiakb\r;
        if(norm(r) < tol)
            fprintf('Converged. Norm(r) : %17.16f\n',norm(r));
            break;
        end
        if(norm(u_Newton) < eps*10)
            fprintf('Newton update is very small: %17.16f\n',...
                norm(u_Newton));
            break;
        end
        if(maxIter < 0)
            fprintf('Maximum iterations exceeded.\n');
            break;
        end
        maxIter = maxIter - 1;
        x(unknownDOFs) = x(unknownDOFs) + u_Newton;        
        %fprintf('Iteration Number: %d norm(r): %17.16f\n',100-maxIter,...
        %    norm(r));
    end    
    X = x;
    H = H.*L;
end

%******************* Uniformity of Deformation Gradient ******************%

% To check for uniform deformation gradient we will caclulate it at
% quadrature point of each element.

F = zeros(3,3,size(IEN,1));
P = zeros(3,3,size(IEN,1));

tmp2 = (1:dofPerNode).'; % We will use it to calculate global DOF number
tmp2 = repmat(tmp2,[size(IEN,2),1]);
Lambda = H./H_orig;
for z=1:size(IEN,1)
    eleNodeNum = IEN(z,:);    
    % Fancy stuff to calculate globalDOF numbers for the global node
    % numbers without using for-loop
    tmp1 = repmat(eleNodeNum,[dofPerNode,1]);
    tmp1 = reshape(tmp1,[numel(tmp1),1]);
    eleGlobalDOF = dofPerNode*(tmp1-1) + tmp2;   
    
    X_ele = X_orig(eleGlobalDOF);
    x_ele = x(eleGlobalDOF);
    [F(:,:,z),P(:,:,z)] = calcDefGrad(x_ele,X_ele,Lambda(z),size(IEN,2),...
        quadOrder,lambda,mu);
end

meanF = mean(F,3);
meanP = mean(P,3);
for z=1:size(IEN,1)
    if(norm(abs(F(:,:,z) - meanF)) > 10^(-12))
        error('Deformation gradient is not uniform!\n');
    end
end

fprintf('mean F11 = %17.16f mean F22 = %17.16f\n',meanF(1,1),meanF(2,2));
fprintf('mean P11 = %17.16f mean P22 = %17.16f\n',meanP(1,1),meanP(2,2));


%***************** Plot the reference and final shapes *******************%

% x = reshape(x,3,[]);
% trisurf(IEN,x(1,:),x(2,:),x(3,:));
% hold on;
% X = reshape(X_orig,3,[]);
% trisurf(IEN,X(1,:),X(2,:),X(3,:)+0.5*ones(1,size(X,2)));
% hold off;
