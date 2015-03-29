function[IEN,Xout]=equiTriMesh(X,level,isPlotOn)
%EQUITRIMESHNODES takes vertices of a triangles as an input and returns a
%connectivity matrix and coordinates of its four "child" triangles a,b,c
%and d as shown below. 'level' controls the number of subdivisions. level 0
%returns just the original triangle. 1 returns the nodes shown below
%
%                          3
%                          /\
%                         /  \
%                        /    \
%                       /      \
%                      /________\
%                     /\        /\
%                    /  \      /  \
%                   /    \    /    \
%                  /      \  /      \
%                 /________\/________\
%                1                    2
% We get the same output as Delaunay triangulation using the cloud of
% points generated in this function without making the connectivity
% ourselves

% Tricky formula to find the number of points

n = 2^(level); % n is the number of parts to divide each side in
numPoints = sum(1:n+1);
Xout = zeros(numPoints,size(X,2));

X1 = X(1,:);
X2 = X(2,:);
X3 = X(3,:);

% Vectors along the sides 1-2 and 1-3 of the triangle respectively
v1 = (X2-X1)/n;
v2 = (X3-X1)/n;

% Calculate the all the nodal points
index = 1;
currMax = n+1;
for j=1:n+1
    start = X1 + (j-1)*v2;
    for i=1:currMax
        Xout(index,:) = start + (i-1)*v1;
        index = index + 1;
    end
    currMax = currMax - 1;
end

% Generating connectivity matrix
numTri = 4^(level);
IEN = zeros(numTri,3);
index = 1;

currIncr = n;
start1 = 1;
end1 = start1 + currIncr;
currIncr = currIncr - 1;
start2 = end1 + 1;
end2 = start2 + currIncr;
currIncr = currIncr - 1;
for i=1:n
    %currNumTri = 2*(end2 - start2 + 1) - 1;
    a = start1;
    b = a + 1;
    c = start2;
    while (b <= end1)
        IEN(index,:) = [a,b,c];
        index = index + 1;
        a = a + 1;
        b = b + 1;
        c = c + 1;
    end
    a = start2;
    b = a + 1;
    c = start1 + 1;
    while (b<=end2)
        IEN(index,:) = [a,c,b];
        index = index + 1;
        a = a + 1;
        b = b + 1;
        c = c + 1;
    end
    start1 = start2;
    end1 = end2;
    start2 = end1 + 1;
    end2 = start2 + currIncr;
    currIncr = currIncr - 1;
end

%IEN = delaunay(Xout(:,1:2));
if(isPlotOn)
    trimesh(IEN,Xout(:,1),Xout(:,2),Xout(:,3));
end

end