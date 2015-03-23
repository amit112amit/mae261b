%% Check Assembly for Consistency

close all; clear; clc;
tic;

%************************** Generate the mesh ****************************%

% For a rectangular mesh
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

% Generate T6 mesh
[IEN,X] = convertToT6Mesh(IEN,X);
numNodes = size(X,1);

trisurf(IEN(:,1:3),X(:,1),X(:,2),X(:,3));
hold on;
scatter3(X(:,1),X(:,2),X(:,3));
hold off;

%****************** Boundary Condition and Meta Arrays *******************%

% Need to know the Global DOF number that has to be constrained.

BC = cell(1,3);
BC{1,1} = [0,0,1; 0,0,1;0,0,1];
BC{1,2} = [1;2;3];
BC{1,3} = [0;0;0];

prescribedDOF = getBCmatrix(BC,X);

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

%******************************** Assembly *******************************%

% External transverse load
f = [0;0;1];

% Quadrature order
quadOrder = 1;

% Thickness of the membrane
H = ones(size(IEN,1),1);

% Elastic constants
lambda = 5*10^8;
mu = 1.5*10^6;

% Set the initial guess for displacement
rng(0);
u = 0.001*rand(size(X));
%u = zeros(size(X)); % Uncomment for zero deformation case

% Apply boundary condition u_z = 0 on the nodes on sides of the square.
rowCol = [ceil(prescribedDOF(:,1)/3),mod(prescribedDOF(:,1),3)];
rowCol(rowCol==0) = 3;
for i=1:size(rowCol,1)
    u(rowCol(i,1),rowCol(i,2)) = prescribedDOF(i,2);
end

X = reshape(X.',[dofPerNode*numNodes,1]);
u = reshape(u.',[dofPerNode*numNodes,1]);

% Thickness stretch
L  = ones(size(IEN,1),1);

x = X + u;

% The exact values
[W_exact,r_exact,K_exact,~] = assemblyT6Quad(X,x,H,f,quadOrder,lambda,...
    mu,IEN,ID,L);

h = norm(x)*logspace(-7,-2,10);
errR = zeros(1,length(h));

% We will use the fact that
for z=1:length(h)
    r_approx = zeros(numel(unknownDOFs),1);
    for i=1:numel(unknownDOFs)
        x(unknownDOFs(i)) = x(unknownDOFs(i)) + h(z);
        [Wp,~,~,~] = assemblyT6Quad(X,x,H,f,quadOrder,lambda,mu,IEN,ID,L);
        x(unknownDOFs(i)) = x(unknownDOFs(i)) - 2*h(z);
        [Wm,~,~,~] = assemblyT6Quad(X,x,H,f,quadOrder,lambda,mu,IEN,ID,L);
        x(unknownDOFs(i)) = x(unknownDOFs(i)) + h(z);
        r_approx(i) = (Wp - Wm)/(2*h(z));
    end
    errR(z) = max(max(abs(r_exact-r_approx)))/norm(r_exact);
    %errR(z) = max(max(abs(r_exact-r_approx)));
end

slope1 = (log(errR(9))-log(errR(7)))/(log(h(9))-log(h(7)));
fprintf('Slope for error in f_internal: %17.16f\n',slope1);

errK = zeros(1,length(h));
for z=1:length(h)
    K_approx = zeros(numel(unknownDOFs),numel(unknownDOFs));
    for m=1:numel(unknownDOFs)
        for n=1:numel(unknownDOFs)      
            x(unknownDOFs(n)) = x(unknownDOFs(n)) + h(z);
            [~,rp,~,~] = assemblyT6Quad(X,x,H,f,quadOrder,...
                lambda,mu,IEN,ID,L);
            x(unknownDOFs(n)) = x(unknownDOFs(n)) - 2*h(z);
            [~,rm,~,~] = assemblyT6Quad(X,x,H,f,quadOrder,...
                lambda,mu,IEN,ID,L);
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
slope2 = (log(errK(9))-log(errK(7)))/(log(h(9))-log(h(7)));
fprintf('Slope for error in K: %17.16f\n',slope2);

toc;
