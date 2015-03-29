%% Consistency test for T6MembraneEle

close all;clear;clc;
tic;
dim = 3;
numNodes = 6;

rng(0);

X = [0.5 0 1 0.25 0.5 0.75;...
     1 0 0 0.433 0 0.433;...
     1 1 1 1 1 1];
x = X + 0.001*rand(3,6);
%x = X;

f = [0;0;1];
quadOrder = 2;
H = 1;
lambda = 5*10^8;
mu = 1.5*10^6;
Lambda = 0.5;

X = reshape(X,[dim*numNodes,1]);
x = reshape(x,[dim*numNodes,1]);

[W,fi,fe,K,~] = T6MembraneEle(X,x,H,f,quadOrder,lambda,mu,Lambda);
fi = reshape(fi,[dim,numNodes]);
%fe = reshape(fe,[dim,numNodes]);
K = reshape(K,[dim,numNodes,dim,numNodes]);

%h = rand(1,1)*10^-5;
h = logspace(-7,-3,15);
errF = zeros(1,length(h));
for z=1:length(h)
    fi_approx = zeros(dim,numNodes);
    for i=1:dim
        for a=1:numNodes
            x(3*(a-1)+i) = x(3*(a-1)+i) + h(z);
            [Wp,~,~,~,~] = T6MembraneEle(X,x,H,f,quadOrder,lambda,mu,...
                Lambda);
            x(3*(a-1)+i) = x(3*(a-1)+i) - 2*h(z);
            [Wm,~,~,~,~] = T6MembraneEle(X,x,H,f,quadOrder,lambda,mu,...
                Lambda);
            x(3*(a-1)+i) = x(3*(a-1)+i) + h(z);

            fi_approx(i,a) = (Wp - Wm)/(2*h(z));
        end
    end
    %errF(z) = max(max(abs(fi-fi_approx)))/norm(fi);
    errF(z) = max(max(abs(fi-fi_approx)));
end
figure(1);
loglog(h,errF);
xlabel('log(h)');
ylabel('log(error)');
slope1 = (log(errF(length(h)-1))-log(errF(length(h)-4)))/...
    (log(h(length(h)-1))-log(h(length(h)-4)));
fprintf('Slope for error in f_internal: %17.16f\n',slope1);

h = logspace(-7,-2,15);
errK = zeros(1,length(h));
for z=1:length(h)
    K_approx = zeros(dim,numNodes,dim,numNodes);
    for i=1:dim
        for a=1:numNodes
            for k=1:dim
                for b=1:numNodes
                    x(3*(b-1)+k) = x(3*(b-1)+k) + h(z);
                    [~,fp,~,~,~] = T6MembraneEle(X,x,H,f,quadOrder,...
                        lambda,mu,Lambda);
                    x(3*(b-1)+k) = x(3*(b-1)+k) - 2*h(z);
                    [~,fm,~,~,~] = T6MembraneEle(X,x,H,f,quadOrder,...
                        lambda,mu,Lambda);
                    x(3*(b-1)+k) = x(3*(b-1)+k) + h(z);
                    fp = reshape(fp,[dim,numNodes]);
                    fm = reshape(fm,[dim,numNodes]);
                    K_approx(i,a,k,b) = (fp(i,a) - fm(i,a))/(2*h(z));
                end
            end
        end
    end
    %errK(z) = max(max(max(max(abs((K - K_approx))))))/...
       max(max(max(max(K))));
    errK(z) = max(max(max(max(abs((K - K_approx))))));
end

% loglog(h,errF,h,errK);
% xlabel('log(h)');
% ylabel('log(error)')
% legend('Log of error in f_internal','Log of error in Kiakb');
% slope2 = (log(errK(length(h)-5))-log(errK(length(h)-10)))/...
%     (log(h(length(h)-5))-log(h(length(h)-10)));
% fprintf('Slope for error in K: %17.16f\n',slope2);

figure(2);
loglog(h,errK);
slope2 = (log(errK(length(h)-1))-log(errK(length(h)-4)))/...
    (log(h(length(h)-1))-log(h(length(h)-4)));
fprintf('Slope for error in K: %17.16f\n',slope2);
toc;