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
colormap(parula(32)); colorbar;
axis equal;
xlim([0.5, 8.5]);
set(gca,'XTick',1:8);
set(gca,'YTick',1:20);
ylabel('FEM Mode #');
xlabel('Test Mode #');
title('Modal Assurance Criterion (MAC) between Test and FEM modes');

%Part b - Calculate the Self and Cross, mass weighted orthogonality.
%How many off diagonals are >0.1? How do they compare to the MAC off
%diagonals? In the x-orth, how many cross terms are >.1? How do those
%compare with the MAC?
%Match modes using cross-orth. Does this pairing agree with MAC? Which
%modes are accurately predicted? (Off diagonal <0.1, diagonal >.9)

SO = corl8(PHItestmn,PHItestmn,M);
figure;
imagesc(SO);
colormap(parula(32)); colorbar;
axis equal;
xlim([0.5, 8.5]);
set(gca,'XTick',1:8);
set(gca,'YTick',1:20);
ylabel('FEM Mode #');
xlabel('Test Mode #');
title('Test Modes Mass-Orthogonality (PHItest*M*PHItest)');

CO = corl8(PHImn,PHItestmn,M);
figure;
imagesc(abs(CO));
colormap(parula(32)); colorbar;
axis equal;
ylim([0.5, 20.5]);
xlim([0.5,8.5]);
set(gca,'YTick',1:8);
set(gca,'XTick',1:20);
ylabel('FEM Mode #');
xlabel('Test Mode #');
title('Absolute Value of Mass weighted Cross-Orthogonality (PHIfem*M*PHItest)');
%Part c - Compute % error in FEM freq, Plot the 8 Test modes with the
%matched FEM modes (same plot). How many modes are accurate (<10%error)?

%Plot FEM and Test mode pairs.
figure;
subplot(4,1,1); hold on;
plot(PHImn(:,1),'r');
plot(PHItestmn(:,1),'k');
ylabel('Amplitude');
xlabel('Mode #');
xlim([1 22]);
grid on;
title('FEM & Test Modes');
legend('FEM Mode #1', 'Test Mode #1');

subplot(4,1,2); hold on;
plot(-1*PHImn(:,3),'r');
plot(PHItestmn(:,2),'k');
ylabel('Amplitude');
xlabel('Mode #');
xlim([1 22]);
grid on;
legend('FEM Mode #3', 'Test Mode #2');
%%%%
% figure;
subplot(4,1,3); hold on;
plot(-1*PHImn(:,2),'r');
plot(PHItestmn(:,3),'k');
ylabel('Amplitude');
xlabel('Mode #');
xlim([1 22]);
grid on;
% title('FEM & Test Modes');
legend('FEM Mode #2', 'Test Mode #3');

subplot(4,1,4); hold on;
plot(PHImn(:,5),'r');
plot(PHItestmn(:,4),'k');
ylabel('Amplitude');
xlabel('Mode #');
xlim([1 22]);
grid on;
legend('FEM Mode #5', 'Test Mode #4');

%%%%
figure;
subplot(4,1,1); hold on;
plot(-1*PHImn(:,6),'r');
plot(PHItestmn(:,5),'k');
ylabel('Amplitude');
xlabel('Mode #');
xlim([1 22]);
grid on;
title('FEM & Test Modes');
legend('FEM Mode #6', 'Test Mode #5');

subplot(4,1,2); hold on;
plot(PHImn(:,8),'r');
plot(PHItestmn(:,6),'k');
ylabel('Amplitude');
xlabel('Mode #');
xlim([1 22]);
grid on;
legend('FEM Mode #8', 'Test Mode #6');
%%%%
% figure;
subplot(4,1,3); hold on;
plot(PHImn(:,11),'r');
plot(PHItestmn(:,7),'k');
ylabel('Amplitude');
xlabel('Mode #');
xlim([1 22]);
grid on;
% title('FEM & Test Modes');
legend('FEM Mode #11', 'Test Mode #7');

subplot(4,1,4); hold on;
plot(PHImn(:,10),'r');
plot(PHItestmn(:,8),'k');
ylabel('Amplitude');
xlabel('Mode #');
xlim([1 22]);
grid on;
legend('FEM Mode #10', 'Test Mode #8');
