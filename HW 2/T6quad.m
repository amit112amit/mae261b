function [N,DN] = T6quad(theta)
%T6QUAD takes parametric coordinates as input and returns the shape
%functions N and shape function derivatives DN for a triangular 6-node
%quadratic element at those coordinates

% The standard domain is the isosceles right angle triangle
%{0<=r<=1,0<=s<=r}. Nodes have been numbered as follows:
% Node 1 -> (1,1) Node 4 -> (0.5,0.5)
% Node 2 -> (0,0) Node 5 -> (0.5,0)
% Node 3 -> (1,0) Node 6 -> (1,0.5)

r = theta(1,1);
s = theta(1,2);

%The Barycentric coordinates
l1 = s;
l2 = 1-r;
l3 = r-s;

%Shape functions
N = [l1*(2*l1-1);l2*(2*l2-1);l3*(2*l3-1);4*l1*l2;4*l2*l3;4*l3*l1];

%Shape function derivatives
DN = [0 4*s-1;...
      4*r-3 0;...
      4*r-4*s-1 1-4*r+4*s;...
      -4*s 4*(1-r);...
      4*(1-2*r+s) 4*(r-1);...
      4*s 4*(r-2*s)];

end