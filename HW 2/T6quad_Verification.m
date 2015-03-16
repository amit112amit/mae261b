%Verification tests for T6quad element

%% Initialization
clear; close all; clc;

tol1 = 10^(-14);

for index=1:1000
    % Create a random point in the domain such that
    % {0<=r<=1,0<=s<=r}
    theta = rand(1,2);
    
    % Get the shape function and derivatives
    [N,DN] = T6quad(theta);
    
    fprintf('Running verification tests for T6quad element...\n');
    %% Partition of unity test for shape functions
    
    if (abs(sum(N) - 1)> tol1)
        error(['Checking Partition of Unity... FAILED! Sum =',...
            '%17.16f'],sum(N));
    else
        fprintf('\tChecking Partition of Unity... PASS\n');
    end
    
    %% Partition of Nullity of shape function derivatives
    
    if (abs(sum(DN(:,1)))>tol1 || abs(sum(DN(:,2)))> tol1)
        error(['Checking Partition of Nullity... FAILED! Sums are ',...
            '%17.16f and %17.16f'],sum(DN(:,1)),sum(DN(:,2)));
    else
        fprintf('\tChecking Partition of Nullity... PASS\n');
    end
    
    %% Consistency Check
    
    tol2 = 10^-6;
    
    h = rand(1)*tol2; %Perturbation in r
    theta(1,1) = theta(1,1) + h;
    [Nplus,~] = T6quad(theta);
    theta(1,1) = theta(1,1) - 2*h;
    [Nminus,~] = T6quad(theta);
    
    dNdr_CentDiff = (Nplus - Nminus)/(2*h);
    theta(1,1) = theta(1,1) + h;
    [~,dNdr] = T6quad(theta);
    
    err1 = abs(dNdr(:,1) - dNdr_CentDiff);
    
    if(err1 > tol2)
        error(['Checking consistency in differentiation with r... FAILED!',...
            ' Error=%17.16f\n'],err1);
    else
        fprintf('\tChecking consistency in differentiation with r... PASS\n');
    end
    
    k = rand(1)*tol2; %Perturbation in s
    theta(1,2) = theta(1,2) + k;
    [Nplus,~] = T6quad(theta);
    theta(1,2) = theta(1,2) - 2*k;
    [Nminus,~] = T6quad(theta);
    
    dNds_CentDiff = (Nplus - Nminus)/(2*k);
    theta(1,2) = theta(1,2) + k;
    [~,dNds] = T6quad(theta);
    
    err2 = abs(dNds(:,2) - dNds_CentDiff);
    
    if(err2 > tol2)
        error(['Checking consistency in differentiation with s... FAILED!',...
            ' Error=%17.16f\n'],err2);
    else
        fprintf('\tChecking consistency in differentiation with s... PASS\n');
    end
    
    %% C0 Completeness Test
    
    tol3 = 10^(-14);
    
    % A linear polynomial in 2-D looks like a*r + b*s + c = 0. The set of
    % coefficients is (a,b,c);
    
    c = rand(1,3);
    p = @(r,s) (c(1)*r + c(2)*s + c(3)); %Anonymous function for polynomial p
    
    theta_star = rand(1,2);
    p_star = p(theta_star(1),theta_star(2));
    
    p_a = [p(1,1);p(0,0);p(1,0);p(0.5,0.5);p(0.5,0);p(1,0.5)];
    [N_a,~] = T6quad(theta_star);
    
    p_Istar = sum(p_a.*N_a);
    
    if(abs(p_Istar - p_star)>tol3)
        error('Checking C0 completeness... FAILED! Error=%17.16f',...
            abs(p_Istar - p_star));
    else
        fprintf('\tChecking C0 completeness... PASS\n');
    end
end
fprintf('Verification of T6quad element completed.\n');