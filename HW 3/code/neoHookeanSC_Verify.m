%% Verification tests for the neoHookeanSC function
clear;
close all;
clc;

F = rand(3);
%We will ensure that F is positive definite so that det(F)>0
F = F*F';
if(det(F)<0)
    error('Invalid deformation gradient generated!');
end
C = F.'*F;
% lambda = rand(1)*10;
% mu = rand(1)*10;
lambda = 5*10^8;
mu = 1.5*10^6;

%Calculate the analytical P and TM
[w_exact,S_exact,TM_exact] = neoHookeanSC(C,lambda,mu);
[w,P,TM] = neoHookean(F,lambda,mu);

errW = abs(w_exact - w)/abs(w);
errorS = max(max(abs(S_exact - F\P)))/norm(F\P);
fprintf('Error in w: %17.16f.\n',errW);
fprintf('Error in S: %17.16f.\n',errorS);


%% Calculate S and TM using 3-point finite difference scheme
h = logspace(-7,-3,100);

errS1 = zeros(1,length(h));
errTM1 = zeros(1,length(h));

for z=1:length(h)
    S_approx = zeros(3,3);
    for i=1:3
        for j=1:3
            C(i,j) = C(i,j) + h(z);
            [Wplus,~,~] = neoHookeanSC(C,lambda,mu);
            C(i,j) = C(i,j)-2*h(z);
            [Wminus,~,~] = neoHookeanSC(C,lambda,mu);
            S_approx(i,j) = (Wplus - Wminus)/h(z);
            C(i,j) = C(i,j) + h(z); %Reset C
        end
    end
    errS1(z)= max(max(abs(S_exact-S_approx)))/norm(S_exact);
    
    %Calculating the tangent modulus numerically
    TM_approx = zeros(3,3,3,3);    
    for I=1:3
        for J=1:3
            for K=1:3
                for L=1:3                    
                    C(K,L) = C(K,L) + h(z);
                    [~,S_plus,~] = neoHookeanSC(C,lambda,mu);
                    C(K,L) = C(K,L) - 2*h(z);
                    [~,S_minus,~] = neoHookeanSC(C,lambda,mu);
                    TM_approx(I,J,K,L) = (S_plus(I,J) - S_minus(I,J))/...
                        (2*h(z));
                    C(K,L) = C(K,L) + h(z);                    
                end
            end
        end
    end    
    errTM1(z) = max(max(max(max((TM_exact - TM_approx)))))/...
        max(max(max(max(TM_exact))));    
end

%% Log-log plot of error vs h
loglog(h,errS1,'r',h,errTM1,'g');
% xlabel('log(h)');
% ylabel('log(error)')
% legend('Log of error in P','Log of error in TM');

slope1 = (log(errS1(length(h)-5))-log(errS1(25)))/(log(h(length(h)-5))-...
    log(h(25)));
fprintf('Slope of log(errS) vs log(h) = %4.3e\n',slope1);
slope2 =(log(errTM1(length(h)-5))-log(errTM1(25)))/(log(h(length(h)-5))-...
    log(h(25)));
fprintf('Slope of log(errTM) vs log(h) = %4.3e\n',slope2);

%% Calculate C_IJKL using equation from lecture notes.

dim =3;

C_IJKL = zeros(dim,dim,dim,dim);
Finv = inv(F);
S = F\P;
for I=1:dim
    for J=1:dim
        for K=1:dim
            for L=1:dim
                for p=1:dim
                    for q=1:dim
                        C_IJKL(I,J,K,L) = C_IJKL(I,J,K,L) +...
                            0.5*Finv(I,p)*Finv(K,q)*(TM(p,J,q,L)-...
                            (p==q)*S(L,J));
                    end
                end
            end
        end
    end
end

errTM2 = max(max(max(max((TM_exact - C_IJKL)))))/...
        max(max(max(max(TM_exact))));
fprintf('Difference in two C_IJKL = %17.16f\n',errTM2);