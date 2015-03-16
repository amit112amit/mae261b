% This script calculates ranks of stiffness matrix for different
% combination of order of quadrature and triangular elements

clear;clc;

lambda = 5*10^8;
mu = 1.5*10^6;

numTrials = 2;

r1 = zeros(1,numTrials);
r2 = zeros(1,numTrials);
r3 = zeros(1,numTrials);
r4 = zeros(1,numTrials);
r5 = zeros(1,numTrials);
r6 = zeros(1,numTrials);


for nT=1:numTrials
    %% Using Triangular 3-node Linear elements
    
    % Generate valid reference and deformed configurations.
    % We check that the surface normal of the triangle points outward so that
    % determinant of F is not negative. We alos ensure that the triangles have
    % good aspect ratio to avoid running into singularities.
    
    dim = 2;
    numNodes = 3;
    
    X = 2*rand(dim,numNodes);
    AR = aspectRatio(X);
    count = 0;
    while(~checkOutwardNormal(X) && ( AR > 2.5))
        count = count + 1;
        %fprintf('Rejected! Count = %d\n',count);
        X = 2*rand(dim,numNodes);
        AR = aspectRatio(X);
    end
    
    %##### Case 1: Zero deformation, first order quadrature
    
    x = X;
    
    [~,~,K1] = T3LinEle(X,x,1,lambda,mu);
    r1(nT) = findStiffnessRank(K1);
    
    %##### Case 2: Random finite deformation, first order quadrature
    
    x = X + 0.1*rand(dim,numNodes);
    ar = aspectRatio(x);
    while(( ar > 2.75))
        x = X + 0.1*rand(dim,numNodes);
        ar = aspectRatio(x);
    end
    
    [~,~,K2] = T3LinEle(X,x,1,lambda,mu);
    r2(nT) = findStiffnessRank(K2);
    
    
    %% Using Triangular 6-node quadratic elements
    
    % Generate valid reference and deformed configurations.
    % We check that the surface normal of the triangle points outward so that
    % determinant of F is not negative. We alos ensure that the triangles have
    % good aspect ratio to avoid running into singularities.
    
    dim = 2;
    numNodes = 6;
    
    X = zeros(dim,numNodes);
    
    X(:,[1,2,3]) = 2*rand(dim,3);
    AR = aspectRatio(X(:,[1,2,3]));
    count = 0;
    while(~checkOutwardNormal(X(:,[1,2,3])) && ( AR > 2.5))
        count = count + 1;
        %fprintf('Rejected! Count = %d\n',count);
        X(:,[1,2,3]) = 2*rand(dim,3);
        AR = aspectRatio(X(:,[1,2,3]));
    end
    %Generating the mid-side nodes to be mid-point of sides
    %In general, the mid-nodes can be anywhere such that a quadratic
    %interpolation is possible with the end-nodes. But we do not want to run
    %lot of tests to ensure correct node-numbering
    X(:,4) = (X(:,1) + X(:,2))/2;
    X(:,5) = (X(:,2) + X(:,3))/2;
    X(:,6) = (X(:,3) + X(:,1))/2;
    
    
    %##### Case 3: Zero deformation, first order quadrature
    
    x = X;
    
    [~,~,K3] = T6QuadEle(X,x,1,lambda,mu);
    r3(nT) = findStiffnessRank(K3);
    
    %##### Case 4: Zero deformation, second-order quadrature
    
    [~,~,K4] = T6QuadEle(X,x,2,lambda,mu);
    r4(nT) = findStiffnessRank(K4);
    
    %##### Case 5: Random finite deformation, first order quadrature
    x = X + 0.01*rand(dim,numNodes);
    
    [~,~,K5] = T6QuadEle(X,x,1,lambda,mu);
    r5(nT) = findStiffnessRank(K5);
    
    %##### Case 6: Random finite deformation, second order quadrature
    
    [~,~,K6] = T6QuadEle(X,x,2,lambda,mu);
    r6(nT) = findStiffnessRank(K6);
end

%% Print a table
fprintf('For zero deformation:\n');
fprintf('\tShape Fn\tQuad Order\tRank\n');
fprintf('\t  Lin\t\t  1\t\t\t  %d\n',(mean(r1)));
fprintf('\t  Quad\t\t  1\t\t\t  %d\n',(mean(r3)));
fprintf('\t  Quad\t\t  2\t\t\t  %d\n',(mean(r4)));

fprintf('\nFor random finite deformation:\n');
fprintf('\tShape Fn\tQuad Order\tRank\n');
fprintf('\t  Lin\t\t  1\t\t\t  %d\n',mean(r2));
fprintf('\t  Quad\t\t  1\t\t\t  %d\n',mean(r5));
fprintf('\t  Quad\t\t  2\t\t\t  %d\n',mean(r6));

