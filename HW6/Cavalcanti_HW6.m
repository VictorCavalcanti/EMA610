clear all; clc; close all;
load('example6.mat');
%thetanom goes k1-k10 ; m1-m6.
%Write sensitivity matrices for K and M as 6x6x16 matrices.
dMim = zeros(6,6,16); %pre-allocate
[PHI0,Lam0,wn0,wnfreq0,~] = getEigSort(Km,Mm);
for i = 11:16
   dMim(i-10,i-10,i) = 1; 
end
dKim = popSensK();


Ltest = (fTest.*(2*pi)).^2; %eigenvalues for test frequencies
% while max(ferror)>=5 || 5
[D,JJ,DD,thetaf]=xupdate(Mm,Km,thetanom,dMim,dKim,Ltest);


M = diag(thetaf(end-5:end));
k = thetaf(1:10);
K = zeros(6);
K = [sum(k(1:4)) -k(2) -k(3) 0 0 -k(4);...
    -k(2) k(2)+k(5) 0 0 0 -k(5);...
    -k(3) 0 k(3)+k(6) 0 0 -k(6);...
    0 0 0 k(7)+k(9) 0 -k(7);...
    0 0 0 0 k(8)+k(10) -k(8);...
    -k(4) -k(5) -k(6) -k(7) -k(8) sum(k(4:8))];
[PHI,Lam,wn,wnfreq,sortindx] = getEigSort(K,M);
ferror = ((wnfreq-fTest)./fTest).*100;
figure;
bar(ferror);
figure;
plot(JJ);