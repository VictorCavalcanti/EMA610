close all; clc; clear all;
load('hw1binfo');

% E = 10^7;
% I = 0.001302;
% p = 0.098;
% l = 40;
% Physical info in metric
A = 0.00016129;
E = 6.894*10^10;
I = 5.4193*10^-10;
p = 2712.63;
l=1.016;
w1 = sqrt(((E*I)/(3*A*p))*((2*pi)/l)^4)/(2*pi);
lam = [ 4.73004074; 7.85320462; 10.9956079; 11.41371655; 17.2787597];
wans = (lam.^2./(2*pi*l^2)).*(E*I/(A*p))^0.5;


phys15 = 0:(2+(2/3)):40;
phys30 = 0:(1+(1/3)):40;
[phi15,v15] = eig(K15,M15);
[phi30,v30] = eig(K30,M30);
wnMAT15 = sqrt(diag(v15));
wnMAT30 = sqrt(diag(v30));
[wnMAT15,w15ind] = sort(wnMAT15,'ascend');
[wnMAT30,w30ind] = sort(wnMAT30,'ascend');
wn15hrz = wnMAT15/(2*pi);
wn30hrz = wnMAT30/(2*pi);

figure; hold on;
% subplot(1,2,1);hold on;
plot(phys15 ,[ 0; phi15(1:2:end,1); 0],'rx-');
plot(phys30,[ 0; phi30(1:2:end,1); 0],'ko-');
title('Mode 1');
ylabel('Amplitude');
xlabel('Physical location of Node (in)');
legend('15 Element Model', '30 Element Model');

figure; hold on;

% subplot(1,2,2); hold on;
plot(phys15 ,[ 0; phi15(1:2:end,2); 0],'rx-');
plot(phys30,[ 0; phi30(1:2:end,2); 0],'ko-');
title('Mode 2');
ylabel('Amplitude');
xlabel('Physical location of Node (in)');
legend('15 Element Model', '30 Element Model');

figure; hold on;
% subplot(1,3,1); hold on;
plot(phys15 ,[ 0; phi15(1:2:end,3); 0],'rx-');
plot(phys30 ,[ 0; phi30(1:2:end,3); 0],'ko-');
title('Mode 3');
ylabel('Amplitude');
xlabel('Physical location of Node (in)');
legend('15 Element Model', '30 Element Model');

figure; hold on;

% subplot(1,3,2); hold on;
plot(phys15 ,[ 0; phi15(1:2:end,4); 0],'rx-');
plot(phys30 ,[ 0; phi30(1:2:end,4); 0],'ko-');
title('Mode 4');
ylabel('Amplitude');
xlabel('Physical location of Node (in)');
legend('15 Element Model', '30 Element Model');

figure; hold on;

% subplot(1,3,3); hold on;
plot(phys15 ,[ 0; phi15(1:2:end,5); 0],'rx-');
plot(phys30, [ 0; phi30(1:2:end,5); 0],'ko-');
title('Mode 5');
ylabel('Amplitude');
xlabel('Physical location of Node (in)');
legend('15 Element Model', '30 Element Model');