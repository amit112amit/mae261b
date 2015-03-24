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
    reshape(mesh_z,numNodes,1)];

%trisurf(IEN,X(:,1),X(:,2),X(:,3));


%****************** Boundary Condition and Meta Arrays *******************%

% Need to know the Global DOF number that has to be constrained. We will
% also constrain the plate to remain stationary in the z-direction

dofPerNode = 3;

% The range of displacement boundary conditions to apply
factor = (-0.25:0.05:4.5).';
u_BC = (factor)*xlim;

constrainedSide1 = X(X(:,2)== 0,:);
constrainedSide2 = X(X(:,1)== xlim,:);
constrainedSide3 = X(X(:,2)== ylim,:);
constrainedSide4 = X(X(:,1)== 0,:);

BC = cell(1,3);
BC{1,1} = [constrainedSide1;constrainedSide2;constrainedSide3;...
    constrainedSide4;X];
BC{1,2} = [2*ones(size(constrainedSide1,1),1);...
    ones(size(constrainedSide2,1),1);2*ones(size(constrainedSide3,1),1);...
    ones(size(constrainedSide4,1),1);3*ones(size(X,1),1)];
BC{1,3} = zeros(size(BC{1,2},1),1); % Dummy Boundary Condition values

prescribedDOF = getBCmatrix(BC,X);

% Generate ID matrix
ID = zeros(numNodes*dofPerNode,2);
ID(:,1) = 1:numNodes*dofPerNode;
unknownDOFs = (1:numNodes*dofPerNode).';
unknownDOFs(prescribedDOF(:,1)) = [];
constrainedDOFs = prescribedDOF(:,1);
globalEqNum = 1;
for i = 1:numel(unknownDOFs)
    ID(unknownDOFs(i),2) = globalEqNum;
    globalEqNum = globalEqNum + 1;
end

% External transverse load
f = [0;0;0];

% Quadrature order
quadOrder = 1;

% Thickness of the membrane
H = 0.01*xlim*ones(size(IEN,1),1);

% The thickness stretch
L = ones(size(IEN,1),1);
L_orig = L;

% Elastic constants
lambda = 5*10^8;
mu = 1.5*10^6;

% Set the initial guess for displacement
rng(0);

u_steps = ceil(abs(u_BC)/0.3) + 5;

% Set tolerance for Newton Iteration
tol = norm(mu)*norm(H)*10^(-12);

meanF11 = zeros(size(u_BC));
meanP11 = zeros(size(u_BC));
meanP22 = zeros(size(u_BC));

X_temp = reshape(X.',[dofPerNode*numNodes,1]);
for q=1:size(u_BC,1)
    % Initial guess for deformed config without BCs is the reference config
    x = X_temp;
    L = L_orig;
    % Actual Boundary Condition values
    prescribedDOF(:,2)= [-1*u_BC(q)*ones(size(constrainedSide1,1),1);...
        u_BC(q)*ones(size(constrainedSide2,1),1);...
        u_BC(q)*ones(size(constrainedSide3,1),1);...
        -1*u_BC(q)*ones(size(constrainedSide4,1),1);zeros(size(X,1),1)];    
    
    %****************************** Assembly *****************************%
    
    u_incr = prescribedDOF(:,2)/u_steps(q);
    
    for i=1:u_steps(q)
        % Apply BCs to the deformed config
        x(constrainedDOFs) = x(constrainedDOFs) + u_incr;
        
        % Newton Iterations
        maxIter = 100;
        while(1)
            [W,r,kiakb,L] = assemblyT3Lin(X_temp,x,H,f,quadOrder,lambda,...
                mu,IEN,ID,L);
            kiakb = sparse(kiakb);
            u_Newton = -kiakb\r;
            if(norm(r) < tol)
                %fprintf('Converged. Norm(r) : %17.16f\n',norm(r));
                break;
            end
            if(norm(u_Newton) < eps*10)
                %fprintf('Newton update is very small: %17.16f\n',...
                %norm(u_Newton));
                break;
            end
            if(maxIter < 0)
                %fprintf('Maximum iterations exceeded.\n');
                break;
            end
            maxIter = maxIter - 1;
            x(unknownDOFs) = x(unknownDOFs) + u_Newton;
        end
        fprintf(['Iterations: %d norm(r): %17.16f norm(u_Newton):',...
            '%17.16f\n'],100-maxIter,norm(r),norm(u_Newton));
    end
    
    %***************** Uniformity of Deformation Gradient ****************%
    
    % To check for uniform deformation gradient we will caclulate it at
    % quadrature point of each element.
    
    F = zeros(3,3,size(IEN,1));
    P = zeros(3,3,size(IEN,1));
    
    tmp2 = (1:dofPerNode).'; % We use it to calculate global DOF number
    tmp2 = repmat(tmp2,[size(IEN,2),1]);
    Lambda = L;
    for z=1:size(IEN,1)
        eleNodeNum = IEN(z,:);
        % Fancy stuff to calculate globalDOF numbers for the global node
        % numbers without using for-loop
        tmp1 = repmat(eleNodeNum,[dofPerNode,1]);
        tmp1 = reshape(tmp1,[numel(tmp1),1]);
        eleGlobalDOF = dofPerNode*(tmp1-1) + tmp2;
        
        temp = X.';
        X_ele = temp(eleGlobalDOF);
        x_ele = x(eleGlobalDOF);
        [F(:,:,z),P(:,:,z)] = calcFandP(x_ele,X_ele,Lambda(z),...
            size(IEN,2),quadOrder,lambda,mu);
    end
    
    meanF = mean(F,3);
    meanP = mean(P,3);
    for z=1:size(IEN,1)
        if(norm(abs(F(:,:,z) - meanF)) > 10^(-12))
            fprintf('Deformation gradient is not uniform %17.16f!\n',...
                norm(abs(F(:,:,z) - meanF)));
        end
    end
    fprintf('q = %d\n',q);
    fprintf('mean F11 = %17.16f mean F22 = %17.16f\n',meanF(1,1),...
        meanF(2,2));
    fprintf('mean P11 = %17.16f mean P22 = %17.16f\n\n',meanP(1,1),...
        meanP(2,2));
    meanF11(q) = meanF(1,1);
    meanP11(q) = meanP(1,1);
    meanP22(q) = meanP(2,2);
    
end

toc;

figure(1);
plot(meanF11(1:76),meanP11(1:76),meanF11(1:76),meanP22(1:76));
xlabel('F11');
ylabel('Component of first Piola-Kirchoff Stress');
title('Stress-strain behaviour for Equi-biaxial strain');
legend('P11','P22','Location','southeast');
%***************** Plot the reference and final shapes *******************%

% x = reshape(x,3,[]);
% trisurf(IEN,x(1,:),x(2,:),x(3,:));
% hold on;
% X = reshape(X_orig,3,[]);
% trisurf(IEN,X(1,:),X(2,:),X(3,:)+0.5*ones(1,size(X,2)));
% hold off;
