% T6QuadEle Verification Script
clear; close all; clc;

tic;

lambda = 5*10^8;
mu = 1.5*10^6;

quadOrder = 1;

%% Generate valid reference and deformed configurations.
% We check that the surface normal of the triangle points outward so that
% determinant of F is not negative. We alos ensure that the triangles have
% good aspect ratio to avoid running into singularities.

dim = 2;
numNodes = 6;

X = zeros(dim,numNodes);

X(:,[1,2,3]) = 2*rand(dim,3);
AR = aspectRatio(X(:,[1,2,3]));
count = 0;
while(~checkOutwardNormal(X(:,[1,2,3])) && ( AR > 2.5))
    count = count + 1;
    fprintf('Rejected! Count = %d\n',count);
    X(:,[1,2,3]) = 2*rand(dim,3);
    AR = aspectRatio(X(:,[1,2,3]));
end
%Generating the mid-side nodes to be mid-point of sides
%In general, the mid-nodes can be anywhere such that a quadratic
%interpolation is possible with the end-nodes. But we do not want to run
%lot of tests to ensure correct node-numbering
X(:,4) = (X(:,1) + X(:,2))/2;
X(:,5) = (X(:,2) + X(:,3))/2;
X(:,6) = (X(:,3) + X(:,1))/2;

x = X + 0.1*rand(dim,numNodes);

% fprintf('Before -> AR: %4.3f\n',AR);
%plotDeformation(X,x);

%% Consistency check
[~,f_exact,K_exact] = T6QuadEle(X,x,quadOrder,lambda,mu);
f_centDiff = zeros(dim,numNodes);

h = norm(x)*logspace(-6,-3);
err_f = zeros(1,length(h));
errK = zeros(1,length(h));

for z=1:length(h)
    for i=1:dim
        for a=1:numNodes
            x(i,a) = x(i,a) + h(z);
            [w_plus,~,~] = T6QuadEle(X,x,quadOrder,lambda,mu);
            x(i,a) = x(i,a) -2*h(z);
            [w_minus,~,~] = T6QuadEle(X,x,quadOrder,lambda,mu);
            f_centDiff(i,a) = (w_plus - w_minus)/(2*h(z));
            x(i,a) = x(i,a) + h(z);
        end
    end
    err_f(z)= max(max(abs(f_exact-f_centDiff)))/norm(f_exact);
    
    % Caculate TM using 3-point finite difference scheme
    K_approx = zeros(dim,numNodes,dim,numNodes);
    index = 1;
    for i=1:dim
        for a=1:numNodes
            for k=1:dim
                for b=1:numNodes
                    x(k,b) = x(k,b) + h(z);
                    [~,fplus,~] = T6QuadEle(X,x,quadOrder,lambda,mu);
                    x(k,b) = x(k,b) - 2*h(z);
                    [~,fminus,~] = T6QuadEle(X,x,quadOrder,lambda,mu);
                    K_approx(i,a,k,b) = (fplus(i,a)-fminus(i,a))/(2*h(z));
                    x(k,b) = x(k,b) + h(z); %Reset x
                end
            end
        end
    end
    errK(z) = max(max(max(max((K_exact - K_approx)))))/...
        max(max(max(max(K_exact))));
end
%% Error plot
loglog(h,err_f,'r',h,errK,'g');
xlabel('log(h)');
ylabel('log(error)')
legend('Log of error in f','Log of error in K');

slope1 = (log(err_f(length(h)-5))-log(err_f(15)))/(log(h(length(h)-5))-...
    log(h(15)));
fprintf('Slope of log(errF) vs log(h) = %4.3e\n',slope1);
slope2 = (log(errK(length(h)-5))-log(errK(15)))/(log(h(length(h)-5))-...
    log(h(15)));
fprintf('Slope of log(errK) vs log(h) = %4.3e\n',slope2);
fprintf('Average error in f:%4.3e\n',mean(err_f));
fprintf('Average error in K:%4.3e\n',mean(errK));

toc;