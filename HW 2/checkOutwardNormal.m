function isConsistent = checkOutwardNormal(X)
%CHECKOUTWARDNORMAL returns true if the triangle formed by the points in X
%has outward surface normal. X should be a 3x3 matrix with each COLUMN
%representing a point. i.e. X = [X1,X2,X3; Y1,Y2,Y3; Z1,Z2,Z3]. If X has
%only 2 rows, the function inserts Z-coordinate to be 0 for all three
%points.

if size(X,1) == 2
    z0 = zeros(1,size(X,2));
    X = [X;z0];
end

A = (X(:,2) - X(:,1))/norm(X(:,2) - X(:,1));
B = (X(:,3) - X(:,1))/norm(X(:,3) - X(:,1));

normal = cross(A,B);
isConsistent = normal > 0;
if(isConsistent)
    isConsistent = 1;
else
    isConsistent = 0;
end
end