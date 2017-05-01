clear all; clc; close all;
load('example6.mat');
%thetanom goes k1-k10 ; m1-m6.
%Write sensitivity matrices for K and M as 6x6x16 matrices.
dMim = zeros(6,6,16); %pre-allocate
for i = 11:16
   dMim(i-10,i-10,i) = 1; 
end
dKim = popSensK();

[D,JJ,DD,thetaf]=xupdate(Mm,Km,thetanom,dMim,dKim);
