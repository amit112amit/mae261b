%% Driver for square membrane deformed under transverse load problem

clear; close all; clc;
tic;

% Flag to control which element to use
useT6QuadEle = false;

%**************************** Generate Mesh ******************************%

% The dimensions for the square membrane in metres.
xlim = 0.1;
ylim = 0.1;
x_incr = 0.005;
y_incr = 0.005;

x_val = 0:x_incr:xlim;
y_val = 0:y_incr:ylim;
z_val = 1;
[mesh_x,mesh_y,mesh_z] = meshgrid(x_val,y_val,z_val);

% Generate the mesh with 3-node triangles
IEN = delaunay(mesh_x,mesh_y);

% Convert the mesh grid points into transpose(x_ia)
numNodes = numel(mesh_x);

% The node co-ordinates for 3-node triangles
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

isScatterOn = false;
plotMesh(IEN(:,1:3),reshape(X.',[],1),'initialMesh.eps',isScatterOn);

%****************** Boundary Condition and Meta Arrays *******************%

% Need to know the Global DOF number that has to be constrained. We will
% also constrain the plate to remain stationary in the z-direction

dofPerNode = 3;

% The Boundary condition cell array
BC = cell(1,3);

% Get co-ordinates of all nodes on the sides of the square membrane
constrainedSide1 = X(X(:,2)== 0,:);
constrainedSide2 = X(X(:,1)== xlim,:);
constrainedSide3 = X(X(:,2)== ylim,:);
constrainedSide4 = X(X(:,1)== 0,:);

temp = [constrainedSide1;constrainedSide2;constrainedSide3;...
    constrainedSide4];

% Remove duplicate rows from the matrix temp
temp = sortrows(temp);
difference = abs(diff([temp;[NaN,NaN,NaN]]));
temp = temp(sum(difference,2)~=0,:);
temp = repmat(temp,[3,1]); % The points to be constrained
BC{1,1} = temp;
BC{1,2} = [ones(size(temp,1)/3,1);2*ones(size(temp,1)/3,1);...
    3*ones(size(temp,1)/3,1)];% The direction to be constrained (z)
BC{1,3} = zeros(size(BC{1,1},1),1); % Constraint value (0)

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
L = ones(size(IEN,1),1);
L_orig = L;

% Elastic constants in N/m^2
lambda = 4*10^6;
mu = 4*10^4;

%****************************** Assembly *****************************%

% Calculate initial guess for displacement in z-direction.
% We are taking advantage of features peculiar to mesh geometry when using
% reshape()
xtemp = reshape(X(:,1),sqrt(numNodes),sqrt(numNodes));
ytemp = reshape(X(:,2),sqrt(numNodes),sqrt(numNodes));
u_z = 10^(-2)*xlim*(sin(xtemp*pi/xlim)+sin(ytemp*pi/ylim));
u_z = reshape(u_z,numNodes,1);

% Set x and y components of displacement to 0
u = [zeros(numNodes,2),u_z];

% Apply boundary condition u_z = 0 on the nodes on sides of the square.
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
    f = repmat([0,0,i*deltaF],size(IEN,1),1);
    
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
fig = figure('visible','off');
plot(f_inc,u_max);
xlabel('Force (N)');
ylabel('Maximum deflection (m)');
title('Force vs. Deflection Curve');
print(fig,'forceDeflection.eps','-depsc');
savefig(fig,'forceDeflection.fig');

toSaveData = [f_inc,u_max];
dlmwrite('forceDeflection.dat',toSaveData,'delimiter','\t',...
    'precision',17);

if(useT6QuadEle)
    isScatterOn = false;
end

isVisible = false;
plotMesh(IEN(:,1:3),x,'deformedMesh.eps',isScatterOn,isVisible);