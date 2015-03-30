%% HW 3 Balloon Problem

close all; clear; clc;
diary('balloonLog.txt');
diary on;
useT6QuadEle = false;

% Generate equilateral triangle mesh
intercept = 10;
X = [intercept,0,0;0,intercept,0;0,0,intercept];
[IEN,X] = equiTriMesh(X,4,false);

% Project the points on plane of equilateral triangle to a sphere
% The radius of the sphere will be same as the intercept value
NormX = sqrt(sum(X.^2,2));
NormX = repmat(NormX,1,3);
X = X./NormX*intercept;

if(useT6QuadEle)
    % Convert the mesh to have 6-node triangles
    [IEN,X] = convertToT6Mesh(IEN,X);
    quadOrder = 2;
else
    quadOrder = 1;
end

%Plot mesh but fileName = empty string, isScatterOn = true, isVisible=true
fig = figure;
fig = plotMesh(fig,IEN,reshape(X.',[],1),'',false,true);

numNodes = size(X,1);
numEle = size(IEN,1);
dofPerNode = 3;

%************************* Boundary Conditions ***************************%

% Get all points in the yz-plane i.e. all points with x coordinate = 0
constrainedEdge1 = X(abs(X(:,1)) < 1e-15,:);

% Get all points in the zx-plane i.e. all points with y coordinate = 0
constrainedEdge2 = X(abs(X(:,2)) < 1e-15,:);

% Get all points in the xy-plane i.e. all points with z coordinate = 0
constrainedEdge3 = X(abs(X(:,3)) < 1e-15,:);

BC{1,1} = [constrainedEdge1; constrainedEdge2; constrainedEdge3];

BC{1,2} = [ones(size(constrainedEdge1,1),1);...
    2*ones(size(constrainedEdge2,1),1);...
    3*ones(size(constrainedEdge3,1),1)];

BC{1,3} = zeros(size(BC{1,2},1),1);

prescribedDOFs = getBCmatrix(BC,X);

% Generate ID matrix
ID = zeros(numNodes*dofPerNode,2);
ID(:,1) = 1:numNodes*dofPerNode;

unknownDOFs = (1:numNodes*dofPerNode).';
unknownDOFs(prescribedDOFs(:,1)) = [];

globalEqNum = 1;
for i = 1:numel(unknownDOFs)
    ID(unknownDOFs(i),2) = globalEqNum;
    globalEqNum = globalEqNum + 1;
end

%************************ Set problem parameters *************************%

% External transverse load in N/m^2
p_max = 1000;
p_min = -20;
p_val = [-1:-1:-20,-15:5:0,10:10:300].';
%deltaP = 10;
%p_steps = p_max/deltaP;
p_steps = numel(p_val);

% Thickness of the membrane in metres
H = 0.001*ones(size(IEN,1),1);

% The thickness stretch
L = 1.2*ones(size(IEN,1),1);

% Elastic constants in N/m^2
lambda = 4*10^6;
mu = 4*10^4;

% Radius in reference configuration
%R_ref = intercept;
R_ref = mean(sqrt(sum(X.^2,2)));

%****************************** Assembly *****************************%

% Set x and y components of displacement to 0
u = 1e-9*rand(size(X));

% Apply boundary condition u_z = 0 on the nodes on sides of the square.
rowCol = [ceil(prescribedDOFs(:,1)/3),mod(prescribedDOFs(:,1),3)];
rowCol(rowCol==0) = 3;
for i=1:size(rowCol,1)
    u(rowCol(i,1),rowCol(i,2)) = prescribedDOFs(i,2);
end

% Calculating pressure along the surface normal of each triangular
% element
X1 = X(IEN(:,1),:); % Coordinates of first vertex
X2 = X(IEN(:,2),:); % Coordinates of second vertex
X3 = X(IEN(:,3),:); % Coordinates of first vertex

% Vectors along two sides of triangle with common vertex 1
side1 = X2 - X1;
side2 = X3 - X1;

% Take cross-product of corresponding rows
crossProd = cross(side1,side2,2);

% Find unit vector of the cross product of two sides to get surface normal
norm_cp = sqrt(sum(crossProd.^2,2));
surfNorm = crossProd./repmat(norm_cp,1,3);

% Prepare column vectors for passing to assembly function
X = reshape(X.',[],1);
u = reshape(u.',[],1);

% Calculate initial guess for deformed configuration
x = X + u;

tol = norm(mu)*norm(H)*10^(-12);
stretchRatio = zeros(p_steps,1);
pressure = zeros(p_steps,1);
for i=1:p_steps
    %fprintf('Pressure step: %d\n',i);    
    % Increment the pressure as force per unit area along the surface
    % normal
    %p = surfNorm*i*deltaP;
    p = surfNorm*p_val(i);
    
    % Tricks to visualize scaled transverse force vectors
%     %quotient = i*deltaP;
%     quotient = p(i);
%     factor = 1;
%     while(quotient >= 10)
%        quotient = quotient/10;
%        factor = factor*10;
%     end
%     scale = 1/factor;
    fig = plotMesh(fig,IEN,x,'',false,true);
%     fig = plotQuiver(fig,p,scale,IEN,x);
    
    % Newton Iterations
    maxIter = 20;
    while(1)
        fprintf('Pressure step: %4.3f Equilibrium Iteration: %d\n',...
            p_val(i),21-maxIter);
        if(useT6QuadEle)
            [W,r,kiakb,L] = assemblyT6Quad(X,x,H,p,quadOrder,lambda,mu,...
                IEN,ID,L);
        else
            [W,r,kiakb,L] = assemblyT3Lin(X,x,H,p,quadOrder,lambda,mu,...
                IEN,ID,L);
        end
        %kiakb = sparse(kiakb);
        u_Newton = -kiakb\r;
        if(norm(r) < tol)
            %fprintf('Converged. Norm(r) : %17.16f\n',norm(r));
            break;
        end
        if(norm(u_Newton) < eps*10)
            %fprintf('Newton update is very small: %17.16f\n',...
            %    norm(u_Newton));
            break;
        end
        if(maxIter < 0)
            %fprintf('Maximum iterations exceeded.\n');
            break;
        end
        maxIter = maxIter - 1;
        x(unknownDOFs) = x(unknownDOFs) + u_Newton;
    end
    if(maxIter <= 0 && norm(r) > 1000*tol)
        diary off;
        t2 = pressure(stretchRatio > 0);
        t1 = stretchRatio(stretchRatio > 0);
        dlmwrite('pressStretch.dat',[t1,t2],'delimiter','\t',...
    'precision',17);
        error('Equilibrium Newton iterations did not converge.');
    end
    fprintf(['Iterations: %d norm(r): %17.16f norm(u_Newton): ',...
        '%17.16f\n'],21-maxIter,norm(r),norm(u_Newton));
    
    % Caculate deformed radius as average of distance of all nodes from the
    % origin. It should be same as one of the intercepts
    temp = reshape(x,3,[]);    
    r_def = mean(sqrt(sum(temp.^2,1)));
    
    stretchRatio(i) = r_def/R_ref; % Ratio of reference and deformed radii
    pressure(i) = p_val(i);
end
toc;

%************************** Plot the results *****************************%
fig = figure('visible','off');
plot(pressure,stretchRatio);
xlabel('Pressure (N/m^2)');
ylabel('Stretch Ratio');
title('Pressure vs. Lateral Stretch Ratio Curve');
print(fig,'pressStretch.eps','-depsc');
savefig(fig,'pressStretch.fig');

toSaveData = [pressure,stretchRatio];
dlmwrite('pressStretch.dat',toSaveData,'delimiter','\t',...
    'precision',17);

if(useT6QuadEle)
    isScatterOn = false;
end

fig = figure;
isVisible = false;
plotMesh(fig,IEN(:,1:3),x,'deformedMeshPr4.eps',isScatterOn,isVisible);
