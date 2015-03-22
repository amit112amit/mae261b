function [N,DN]=T3Lin(theta)
%T3LIN takes parametric coordinates as input and returns the shape
%functions N and shape function derivatives DN for a triangular three node
%linear element at those coordinates. 'theta' must be 1-by-2 - first
%element is r coordinate and second is s coordinate in the standard
%domain

% The standard domain is the isosceles right angle triangle
%{0<=r<=1,0<=s<=r}. Nodes have been numbered as follows:
% Node 1 -> (1,1)
% Node 2 -> (0,0)
% Node 3 -> (1,0)

r = theta(1,1);
s = theta(1,2);

%The Barycentric coordinates
l1 = s;
l2 = 1-r;
l3 = r-s;

%Shape functions
N = [l1;...
     l2;...
     l3];

%Shape function derivatives
DN = [0,1;...
      -1,0;...
      1,-1];

end