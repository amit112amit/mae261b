% In this script we will use plane stress neo-Hookean implementation to
% solve for stresses for given uni-axial and equi-bi-axial deformation

clear;close all;clc;

% Stretch ratio should be greater than 0. It should be less than 1 for
% compression and greater than 1 for elongation
F11 = linspace(0.4,10,1000);
%F11 = [F11,linspace(1,2.5,100)];
%F11 = [F11,linspace(2.5,5,50)];

%We will use properties of rubber
lambda = 5*10^8;
mu = 1.5*10^6;

%% Uniaxial deformation
F_uni = zeros(2,2,length(F11));
P_uni = zeros(2,2,length(F11));
guessF33 = ones(1,length(F11));
for i=1:length(F11)
    F_uni(:,:,i) = [F11(i), 0; 0, 1];
    temp = squeeze(F_uni(:,:,i));
    [~,P_uni(:,:,i),~,F33] = planeStressNH(temp,lambda,mu,guessF33(i));
    if(i~=1000)
        guessF33(i+1)=F33;
    end
end

% Plotting stress-strain curve.
temp1 = squeeze(P_uni(1,1,:));
temp2 = squeeze(P_uni(2,2,:));
figure(1);
plot(F11,temp1,'r',F11,temp2,'b','LineWidth',2);
xlabel('F11');
ylabel('Component of Piola-Kirchoff Stress');
legend('P11','P22','Location','southeast');
title('Stress-strain Behavior under Uniaxial Strain');
figure(3);
plot(F11,guessF33);

%% Equibiaxial deformation
F_bi = zeros(2,2,length(F11));
P_bi = zeros(2,2,length(F11));
prevF33 = 1;
for i=1:length(F11)
    F_bi(:,:,i) = [F11(i), 0; 0, F11(i)];
    temp = squeeze(F_bi(:,:,i));
    [~,P_bi(:,:,i),~,F33] = planeStressNH(temp,lambda,mu,prevF33);
    prevF33 = F33;
end

% Plotting stress-strain curve.
temp1 = squeeze(P_bi(1,1,:));
temp2 = squeeze(P_bi(2,2,:));
figure(2);
plot(F11,temp1,'--r',F11,temp2,'-.b','LineWidth',2);
xlabel('F11');
ylabel('Component of Piola-Kirchoff Stress');
legend('P11','P22','Location','southeast');
title('Stress-strain Behavior under Equibiaxial Strain');
