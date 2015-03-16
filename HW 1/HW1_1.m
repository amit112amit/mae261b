%************************ MAE 261B HW 1 Part I **********************%
%Curvilinear Kinematics:
%This script calculates 
%1. Tangent and dual basis vectors in reference and deformed 
%   configurations
%2. Deformation gradient, F
%3. Cauchy-Green deformation tensor, C
%4. Green Strain, E
%
%Author: Amit Singh
%Year: 2014
%********************************************************************%
clear;
close all;
clc;

% Initialize parameters
a0 = 1;
lambda1 = 1;
lambda2 = 2;

R = a0;
phi = pi/4;

% The cylindrical basis vectors
e_R = [cos(phi), sin(phi), 0]';
e_phi = [-sin(phi), cos(phi), 0]';
e_Z = [0,0,1]';

% The basis vectors in reference configurations - Gt for tangent basis
% vectors and Gd for dual basis vectors
Gt =[e_R, R*e_phi, e_Z];
Gd = [e_R, e_phi/R,e_Z];

% The basis vectors in reference configurations - gt for tangent basis
% vectors and gd for dual basis vectors
gt = [lambda1*e_R, lambda1*R*e_phi, lambda2*e_Z];
gd = [(1/lambda1)*e_R, (1/(lambda1*R))*e_phi, (1/lambda2)*e_Z];

%Deformation Gradient, $ F = gt(:,i) \otimes Gd(:,i) $
F = gt(:,1)*Gd(:,1)' + gt(:,2)*Gd(:,2)' + gt(:,3)*Gd(:,3)';

%The right Cauchy-Green deformation tensor
C = F'*F;

%The Green strain
E = (C - eye(3))/2;

display(gt);
display(gd);
display(Gt);
display(Gd);
display(F);
display(C);
display(E);
