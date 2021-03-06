%% Testing T6 mesh generation

close all; clear; clc;
useT6QuadEle = true;
%************************** Generate the mesh ****************************%

% For rectangular mesh
xlim = 1;
ylim = 1;
zlim = 1;
[mesh_x,mesh_y,mesh_z] = meshgrid(0:xlim,0:ylim,1:zlim);

% Generate the mesh
IEN = delaunay(mesh_x,mesh_y);

% Convert the mesh grid points into x_ia
numNodes = numel(mesh_x);

% The reference configuration
X = [reshape(mesh_x,numNodes,1), reshape(mesh_y,numNodes,1),...
    reshape(mesh_z,numNodes,1)];

if(useT6QuadEle)
    % Convert the mesh to have 6-node triangles
    [IEN,X] = convertToT6Mesh(IEN,X);
    numNodes = size(X,1);
    quadOrder = 2;
else
    quadOrder = 1;
end

trisurf(IEN(:,1:3),X(:,1),X(:,2),X(:,3));
hold on;
scatter3(X(:,1),X(:,2),X(:,3));
hold off;

%****************** Boundary Condition and Meta Arrays *******************%

% Need to know the Global DOF number that has to be constrained. We will
% also constrain the plate to remain stationary in the z-direction

dofPerNode = 3;

% The Boundary condition cell array
BC = cell(1,3);

BC{1,1} = [0,0,1;0,0,1;0,0,1]; % The points to be constrained
BC{1,2} = [1;2;3]; % The direction to be constrained (z)
BC{1,3} = [0.001;0.001;0.001]; % Constraint value (0)

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

%************************ Set problem parameters *************************%

% External transverse load in N/m^2
f_max = 1000;
deltaF = 10;
f_steps = f_max/deltaF;

% Thickness of the membrane in metres
H = 0.001*ones(size(IEN,1),1);

% The thickness stretch
L = 0.5*ones(size(IEN,1),1);

% Elastic constants in N/m^2
lambda = 4*10^6;
mu = 4*10^4;

%****************************** Assembly *****************************%

% Back-up original configuration
H_orig = H;
X_orig = X;

u_z = 0.01*rand(numNodes,1);

% Set x and y components of displacement to 0
u = [zeros(numNodes,2),u_z];

% Apply boundary condition u_z = 0 on constrained DOF
rowCol = [ceil(prescribedDOF(:,1)/3),mod(prescribedDOF(:,1),3)];
rowCol(rowCol==0) = 3;
for i=1:size(rowCol,1)
    u(rowCol(i,1),rowCol(i,2)) = prescribedDOF(i,2);
end

% Prepare for passing to assembly function
X = reshape(X.',[dofPerNode*numNodes,1]);
u = reshape(u.',[dofPerNode*numNodes,1]);

% Calculate initial guess for deformed configuration
x = X + u;

tol = norm(mu)*norm(H)*10^(-12);
u_max = zeros(f_steps,1);
f_inc = zeros(f_steps,1);
for i=1:f_steps
    % Increment the force
    f = [0;0;i*deltaF];
    
    % Newton Iterations
    maxIter = 100;
    while(1)
        if(useT6QuadEle)
            [W,r,kiakb,L] = assemblyT6Quad(X,x,H,f,quadOrder,lambda,mu,...
                IEN,ID,L);
        else
            [W,r,kiakb,L] = assemblyT3Lin(X,x,H,f,quadOrder,lambda,mu,...
                IEN,ID,L);
        end
        %kiakb = sparse(kiakb);
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
    end
    fprintf(['Iterations: %d norm(r): %17.16f norm(u_Newton): ',...
        '%17.16f\n'],100-maxIter,norm(r),norm(u_Newton));
    
    % Caculate displacement
    u = x - X;
    u_max(i) = max(abs(u)); % Maximum deflection
    f_inc(i) = norm(f);
end
toc;

%************************** Plot the results *****************************%
figure(1)
plot(f_inc,u_max);
xlabel('Force (N)');
ylabel('Maximum deflection (m)');
title('Force vs. Deflection Curve');

figure(2);
x = reshape(x,[dim,numNodes]);
surf(x(1,:),x(2,:),x(3,:));