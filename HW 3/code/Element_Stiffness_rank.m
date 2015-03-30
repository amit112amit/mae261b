% This script calculates ranks of stiffness matrix for different
% combination of order of quadrature and triangular elements

clear;clc;
f = [0;0;1];
H = 1;
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
    dim = 3;
    numNodes = 3;
    X = [0.5 0 1; 1 0 0; 1 1 1];   
    
    %##### Case 1: Zero deformation, first order quadrature    
    x = X;    
    [~,~,~,K1] = T3MembraneEle(X,x,H,f,1,lambda,mu,1);    
    r1(nT) = findStiffnessRank(K1);
    
    %##### Case 2: Random finite deformation, first order quadrature    
    x = X + 0.001*rand(dim,numNodes);    
    [~,~,~,K2] = T3MembraneEle(X,x,H,f,1,lambda,mu,1);
    r2(nT) = findStiffnessRank(K2);    
    
    %% Using Triangular 6-node quadratic elements    
    dim = 3;
    numNodes = 6;    
    X = [0.5 0 1 0.25 0.5 0.75;...
        1 0 0 0.433 0 0.433;...
        1 1 1 1 1 1];
        
    %##### Case 3: Zero deformation, first order quadrature    
    x = X;    
    [~,~,~,K3] = T6MembraneEle(X,x,H,f,1,lambda,mu,1);
    r3(nT) = findStiffnessRank(K3);
    
    %##### Case 4: Zero deformation, second-order quadrature    
    [~,~,~,K4] = T6MembraneEle(X,x,H,f,2,lambda,mu,1);
    r4(nT) = findStiffnessRank(K4);
    
    %##### Case 5: Random finite deformation, first order quadrature
    x = X + 0.001*rand(dim,numNodes);    
    [~,~,~,K5] = T6MembraneEle(X,x,H,f,1,lambda,mu,1);
    r5(nT) = findStiffnessRank(K5);
    
    %##### Case 6: Random finite deformation, second order quadrature    
    [~,~,~,K6] = T6MembraneEle(X,x,H,f,2,lambda,mu,1);
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

