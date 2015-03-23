function [BC]  = getBCmatrix(X_BC,X)
%GETBCMATRIX takes the boundary condition in form of the co-ordinates of
%the node the DOF being constrained and the value of prescribed boundary
%condition. X_BC is always a cell array of size 1x3.
% Eg: X_BC = {[2,2,1],[3],[0]} means that the node at the origin i.e.
% [0,0,0] is constrained along the z-axis (indicated by 3) has 0 
% displacement. X contains coordinates of all the nodes in the system
% arranged in order of their global node numbers. The structure of X
% should be like
%
% X = transpose([1 1 2 0 3 1 4;...
%      0 2 2 0 1 2 2;...
%      1 1 1 1 1 1 1]);
% Here, there are 7 nodes and coordinates of third node are (2,2,1)

dofPerNode = 3;

nodeCoordinates = X_BC{1,1};
nodeDOF = X_BC{1,2};
nodeBC = X_BC{1,3};

BC = zeros(size(nodeCoordinates,1),2);

for i=1:size(nodeCoordinates,1)
    currNode = nodeCoordinates(i,:);
    currNodeLocalDOF = nodeDOF(i);
    currBC = nodeBC(i);
    
    currNode = repmat(currNode,[size(X,1),1]);    
    globalNodeNum = currNode - X;
    exactMatch = globalNodeNum == 0;
    exactMatch = sum(exactMatch,2);
    [globalNodeNum,~] = find(exactMatch == 3);
    
    globalDOFNum = dofPerNode*(globalNodeNum - 1) + currNodeLocalDOF;
    BC(i,:) = [globalDOFNum,currBC];
end