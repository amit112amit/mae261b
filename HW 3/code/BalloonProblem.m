%% HW 3 Balloon Problem

close all; clear; clc;

useT6QuadEle = true;

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
    [IEN,X] = convertToT6Mesh(IEN,X);
end

%Plot mesh but fileName = empty string, isScatterOn = true, isVisible=true
plotMesh(IEN,reshape(X.',[],1),'',true,true);

%************************* Boundary Conditions ***************************%