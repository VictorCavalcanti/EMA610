clear all; clc; close all;
load('beam.mat');
%COMPUTE MASS NORMALIZED MODES FROM beam.mat
[PHI,LAM] = eig(K,M); %eig problem
wnhrz = sqrt(abs(diag(LAM)))./(2*pi); %freq in hrtz.
[wnhrz, ascindx] = sort(wnhrz,'ascend'); %Sort frequencies, ascend.
PHI = PHI(:,ascindx); %sort shapes based on freqs.
mnorm = zeros(22,1);
PHImn = zeros(22,22);
for i = 1:22
 mnorm(i) = 1/(sqrt(PHI(:,i)'*M*PHI(:,i)));
 PHImn(:,i) = mnorm(i)*PHI(:,i);
end
%REMOVE RIGID BODY MODES
PHImn = PHImn(:,3:end);

%IMPORT DATA FROM 'test.dat'

A = importdata('test.dat');
wtest = A(1:8);
PHItest = zeros(22,8);
for i = 1:22
PHItest(i,1:6) = A(9+6*i-5:9+6*i);
PHItest(i,7:8) = A(142+2*i-1:142+2*i);
end
testnorm = zeros(8,1);
PHItestmn = zeros(22,8);
for i = 1:8
 testnorm(i) = 1/(sqrt(PHItest(:,i)'*M*PHItest(:,i)));
 PHItestmn(:,i) = testnorm(i)*PHItest(:,i);
end

%Part a - Calculate the MAC matrix between the test modes and the FEM
%modes
%Mode matching pairs (TEST, FEM): (1,1) (2,3) (3,2) (4,5) (5,6) (6,8)
%(7,11) (8,10)
MAC = mac(PHImn,PHItestmn); %Call MAC function
%Color plot MAC.
imagesc(MAC);
colormap(jet(100)); colorbar;
axis equal;
xlim([0.5, 8.5]);
set(gca,'XTick',1:8);
set(gca,'YTick',1:20);
ylabel('FEM Mode #');
xlabel('Test Mode #');
title('Modal Assurance Criterion (MAC) between Test and FEM modes');

%Part b - Calculate the Self and Cross, mass weighted orthogonality.
SO = corl8(PHItestmn,PHItestmn,M);
CO = corl8(PHItestmn,PHImn,M);










