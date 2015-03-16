% Verification of TriGaussQuad function

clear; close; clc;

fprintf('Running verification tests for Gauss Quadrature...\n');

tol = 10^-15;

for z=1:1000 % Regression testing
    %% First order rule.

    % A random linear polynomial in (r,s) has the form (a*r + b*s + c)
    % Integrating (a*r+b*s+c) over the domain {0<=r<=1,0<=s<=r} gives
    % (a/3 + b/6 + c/2)
    coeff = rand(1,3);
    p1 = @(r,s) (coeff(1)*r + coeff(2)*s + coeff(3));
    
    exact_integ = coeff(1)/3 + coeff(2)/6 + coeff(3)/2;
    
    [points,weights] = TriGaussQuad(1);
    
    gauss_integ = 0;
    for index = 1:length(weights)
        gauss_integ = gauss_integ + ...
            p1(points(index,1),points(index,2))*weights(index);
    end
    
    err1 = abs(exact_integ - gauss_integ);
    
    if( err1 > tol)
        error(['Checking linear triangular Gauss Quadrature... FAILED!',...
            ' Error = %17.16f'],err1);
    else
        fprintf('\tChecking linear triangular Gauss Quadrature... PASS\n');
    end
    
    %% Second Order rule   
    
    % A general quadratic polynomial in (r,s) has the form
    % (a*r^2 + b*s^2 + c*r*s + d*r +e*s + f)
    % Integrating this over the domain {0<=r<=1,0<=s<=r} gives
    % (1/24)*(6*a + 2*b + 3*c + 8*d + 4*e + 12*f)
    coeff = rand(1,6);
    p2 = @(r,s) (coeff(1)*r^2 + coeff(2)*s^2 + coeff(3)*r*s +...
        coeff(4)*r + coeff(5)*s +coeff(6));
    
    exact_integ = (1/24)*(6*coeff(1) + 2*coeff(2) + 3*coeff(3) +...
        8*coeff(4) + 4*coeff(5) + 12*coeff(6));
    
    [points,weights] = TriGaussQuad(2);
    
    gauss_integ = 0;
    for index = 1:length(weights)
        gauss_integ = gauss_integ + ...
            p2(points(index,1),points(index,2))*weights(index);
    end
    
    err2 = abs(exact_integ - gauss_integ);
    
    if( err2 > tol)
        error(['Checking quadratic triangular Gauss Quadrature... FAILED!',...
            ' Error = %17.16f'],err2);
    else
        fprintf('\tChecking quadratic triangular Gauss Quadrature... PASS\n');
    end
end

fprintf('Verification of triangular Gauss Quadrature completed.\n');