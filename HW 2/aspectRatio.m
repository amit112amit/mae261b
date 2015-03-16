function AR = aspectRatio(X)
%ASPECTRATIO returns the aspect ratio of triangle formed by points in X.
%X should have three columns where each column gives a point of the
%triangle. So x-coordinates of all three points must be in same row.
%Similarly for y and z. Formula used to calculate aspect ratio compares
%circumradius and in radius of the triangle.

% Side lengths
a = norm(X(:,2) -X(:,1));
b = norm(X(:,3) -X(:,1));
c = norm(X(:,2) -X(:,3));

s=(a+b+c)/2;
AR = (a*b*c)/(8*(s-a)*(s-b)*(s-c));

end