%% Consistency check for plane stress
clear;
close all;
clc;

tic;

F = rand(2);
%We will ensure that F is positive definite so that det(F)>0
F = F*F';
if(det(F)<0)
    error('Invalid deformation gradient generated!');
end
% lambda = rand(1)*10;
% mu = rand(1)*10;
lambda = 5*10^8;
mu = 1.5*10^6;

guessF33 = 1;
%Calculate the analytical P and TM
[w_exact,P_exact,TM_exact,F33] = planeStressNH(F,lambda,mu,guessF33);

%% Calculate P and TM using 3-point finite difference scheme
h = norm(F)*logspace(-6,-3);

errP = zeros(1,length(h));
errTM = zeros(1,length(h));

for z=1:length(h)
    P_approx = zeros(2,2);
    for i=1:2
        for j=1:2
            F(i,j) = F(i,j) + h(z);
            [Wplus,~,~,~] = planeStressNH(F,lambda,mu,guessF33);
            F(i,j) = F(i,j)-2*h(z);
            [Wminus,~,~,~] = planeStressNH(F,lambda,mu,guessF33);
            P_approx(i,j) = (Wplus - Wminus)/(2*h(z));
            F(i,j) = F(i,j) + h(z); %Reset F
        end
    end
    errP(z)= max(max(abs(P_exact-P_approx)))/norm(P_exact);
    
    % Caculate TM using 3-point finite difference scheme
    TM_approx = zeros(2,2,2,2);
    %index = 1;
    for i=1:2
        for J=1:2
            for k=1:2
                for L=1:2
                    F(k,L) = F(k,L) + h(z);
                    [~,Pplus,~,~] = planeStressNH(F,lambda,mu,guessF33);
                    F(k,L) = F(k,L) - 2*h(z);
                    [~,Pminus,~,~] = planeStressNH(F,lambda,mu,guessF33);
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


%% Error plot
loglog(h,errP,'r',h,errTM,'g');
xlabel('log(h)');
ylabel('log(error)')
legend('Log of error in P','Log of error in TM');

slope1 = (log(errP(length(h)-5))-log(errP(15)))/(log(h(length(h)-5))-...
    log(h(15)));
fprintf('Slope of log(errP) vs log(h) = %4.3e\n',slope1);
slope2 = (log(errTM(length(h)-5))-log(errTM(15)))/(log(h(length(h)-5))-...
    log(h(15)));
fprintf('Slope of log(errTM) vs log(h) = %4.3e\n',slope2);

fprintf('Average error in P:%4.3e\n',mean(errP));
fprintf('Average error in TM:%4.3e\n',mean(errTM));

toc;