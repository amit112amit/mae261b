function [points,weights]= TriGaussQuad(order)
% TRIGAUSSQUAD function returns the sampling points and weights for
% Gauss-Legendre quadrature over a triangular domain. Input order can be
% 1 or 2 for linear and quadratic quadrature respectively. The output
% points is in terms of (r,s) co-ordinates and not the barycentric
% coordiantes

switch order
    case 1
        points = [2/3,1/3];
        weights = 0.5;
    case 2
        points = [0.5 0.5; 1 0.5; 0.5 0];
        weights = [1/6; 1/6; 1/6];
    otherwise
        error('Gauss quadrature for order %d has not been implemented.',...
            order);
end

end