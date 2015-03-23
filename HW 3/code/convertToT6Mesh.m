function [IEN_T6,X_T6] = convertToT6Mesh(IEN,X)
%CONVERTTOT6MESH Converts T3 mesh to T6 mesh
%   Adds mid-side nodes to the triangles

dim = 3;
nodesPerEle = 6;
numEle = size(IEN,1);

X1 = X(IEN(:,1),:);
X2 = X(IEN(:,2),:);
X3 = X(IEN(:,3),:);

X4 = (X1+X2)/2;
X5 = (X2+X3)/2;
X6 = (X3+X1)/2;

X_T6 = [X;X4;X5;X6];

% Remove duplicate rows from X_T6
X_T6 = sortrows(X_T6);
difference = abs(diff([X_T6;[NaN,NaN,NaN]]));
X_T6 = X_T6(sum(difference,2)~=0,:);

X_temp = zeros(numEle,dim,nodesPerEle);
X_temp(:,:,1) = X1;
X_temp(:,:,2) = X2;
X_temp(:,:,3) = X3;
X_temp(:,:,4) = X4;
X_temp(:,:,5) = X5;
X_temp(:,:,6) = X6;

IEN_T6 = zeros(numEle,nodesPerEle);
for i=1:numEle
    for j=1:nodesPerEle
        currX = X_temp(i,:,j);
        currX = repmat(currX,[size(X_T6,1),1]);
        currX = sum(abs(X_T6 - currX),2);
        [globalNodeNum,~] = find(currX == 0);
        IEN_T6(i,j) = globalNodeNum;
    end
end

end

