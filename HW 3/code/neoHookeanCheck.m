%% Verification tests for the neoHookean function
clear;
close all;
clc;

F = rand(3);
%We will ensure that F is positive definite so that det(F)>0
F = F*F';
if(det(F)<0)
   error('Invalid deformation gradient generated!');
end
lambda = rand(1)*10;
mu = rand(1)*10;


%Calculate the analytical P and TM
[w_exact,P_exact,TM_exact] = neoHookean(F,lambda,mu);

%% Calculate P and TM using 3-point finite difference scheme
h = norm(F)*logspace(-6,-3);

errP = zeros(1,length(h));
errTM = zeros(1,length(h));

for z=1:length(h)
    P_approx = zeros(3,3);
    for i=1:3
        for j=1:3
            F(i,j) = F(i,j) + h(z);
            [Wplus,~,~] = neoHookean(F,lambda,mu);
            F(i,j) = F(i,j)-2*h(z);
            [Wminus,~,~] = neoHookean(F,lambda,mu);
            P_approx(i,j) = (Wplus - Wminus)/(2*h(z));
            F(i,j) = F(i,j) + h(z); %Reset F
        end
    end
    errP(z)= max(max(abs(P_exact-P_approx)))/norm(P_exact);
    
    % Caculate TM using 3-point finite difference scheme
    TM_approx = zeros(3,3,3,3);
    index = 1;
    for i=1:3
        for J=1:3
            for k=1:3
                for L=1:3
                    F(k,L) = F(k,L) + h(z);
                    [~,Pplus,~] = neoHookean(F,lambda,mu);
                    F(k,L) = F(k,L) - 2*h(z);
                    [~,Pminus,~] = neoHookean(F,lambda,mu);
                    TM_approx(i,J,k,L) = (Pplus(i,J)-Pminus(i,J))/(2*h(z));
                    F(k,L) = F(k,L) + h(z); %Reset F
                end
            end
        end
    end
    errTM(z) = max(max(max(max((TM_exact - TM_approx)))))/...
        max(max(max(max(TM_exact))));
    %     fprintf('Average error in elements of TM:%4.3e\n',err);
end

%% Log-log plot of error vs h
loglog(h,errP,'r',h,errTM,'g');
xlabel('log(h)');
ylabel('log(error)')
legend('Log of error in P','Log of error in TM');

slope1 = (log(errP(length(h)-5))-log(errP(5)))/(log(h(length(h)-5))-...
    log(h(5)));
fprintf('Slope of log(errP) vs log(h) = %4.3e\n',slope1);
slope2 = (log(errTM(length(h)-5))-log(errTM(5)))/(log(h(length(h)-5))-...
    log(h(5)));
fprintf('Slope of log(errTM) vs log(h) = %4.3e\n',slope2);

%% Checking Material Frame Indifference

v = rand(3,1);
n = v/norm(v);
theta = rand(1)*pi;
n_hat = [0,-n(3),n(2); n(3),0,-n(1);-n(2),n(1),0];
Q = eye(3) + sin(theta)*n_hat + (1- cos(theta))*(n*n' - eye(3));
tol = 10^-7;

% Check if Q belongs to SO3
if (abs(det(Q)-1) > tol || norm(abs(inv(Q) - Q')) > tol)
    error('Q does not belong to SO3.\n');
end

F_rot = Q*F;

[w_rot,P_rot,TM_rot] = neoHookean(F_rot,lambda,mu);

fprintf('Relative error in w due to coordinate rotation = %4.3e\n',...
    abs(w_rot - w_exact)/norm(w_exact));
fprintf('Relative error in P due to coordinate rotation = %4.3e\n',...
    max(max((abs(P_rot - Q*P_exact))))/norm(Q*P_exact));
TM_rot_exact = zeros(3,3,3,3);
for i=1:3
    for J=1:3
        for k=1:3
            for L=1:3
                for m=1:3
                    for n=1:3
                        TM_rot_exact(i,J,k,L) = ...
                            TM_rot_exact(i,J,k,L) + Q(i,m)*Q(k,n)*...
                            TM_exact(m,J,n,L);
                    end
                end
            end
        end
    end
end
fprintf('Relative error in TM due to coordinate rotation = %4.3e\n',...
    max(max(max(max(abs(TM_rot - TM_rot_exact)))))/...
    max(max(max(max(TM_rot_exact)))));

%% Symmetry test
F_sym = F*Q;

[w_sym,P_sym,TM_sym] = neoHookean(F_sym,lambda,mu);

fprintf('Relative error in w for symmetry test = %4.3e\n',...
    abs(w_sym - w_exact)/norm(w_exact));
fprintf('Relative error in P for symmetry test = %4.3e\n',...
    max(max((abs(P_sym - P_exact*Q))))/norm(P_exact*Q));

TM_sym_exact = zeros(3,3,3,3);
for i=1:3
    for J=1:3
        for k=1:3
            for L=1:3
                for M=1:3
                    for N=1:3
                        TM_sym_exact(i,J,k,L) = ...
                            TM_sym_exact(i,J,k,L) + Q(M,J)*Q(N,L)*...
                            TM_exact(i,M,k,N);
                    end
                end
            end
        end
    end
end

fprintf('Relative error in TM for symmetry test = %4.3e\n',...
    max(max(max(max(abs(TM_sym - TM_sym_exact)))))/...
    max(max(max(max(TM_sym_exact)))));

