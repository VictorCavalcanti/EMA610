clear all; clc; close all;
load('beam.mat');

%%
%1a - Compute Mode Shapes and frequencies in Hertz, sort both in ascending,
%plot freq(hrz) vs. mode #. Plot first two elastic modes vs. mode #.
%PHI - Normalized to max unit deflection phi = phi/max(phi).
[PHI,LAM] = eig(K,M); %eig problem
wnhrz = sqrt(abs(diag(LAM)))./(2*pi); %freq in hrtz.
[wnhrz, ascindx] = sort(wnhrz,'ascend'); %Sort frequencies, ascend.
PHI = PHI(:,ascindx); %sort shapes based on freqs.
plot(1:numel(diag(LAM)),wnhrz,'o-');
ylabel('Frequency (Hrz)');
xlabel('Mode #');
title('Frequency (Hrz) vs Mode #');
xlim([1 22]);
grid on;

figure; 
subplot(2,1,1); hold on;
plot(PHI(:,3),'r-');
ylabel('Amplitude');
xlabel('Mode #');
xlim([1 22]);
grid on;
title('Mode Shape - 1st Elastic Mode');
subplot(2,1,2); hold on;
plot(PHI(:,4),'b-');
ylabel('Amplitude');
xlabel('Mode #');
xlim([1 22]);
grid on;
title('Mode Shape - 2nd Elastic Mode');
%%
%1b - It appears that Matlab normalized these modes to maximum unit
%deflection, since the mode shapes maxima are at +-1.
%%
%1c - Normalize the modes to unit length (phi^t*phi=1). Plot phi^t*M*phi &
%phi^t*K*phi vs freq (hrz).
ulnorm = zeros(22,1);
PHIul = zeros(22,22);
for i = 1:22
    ulnorm(i) = 1/(sqrt(PHI(:,i)'*PHI(:,i)));
    PHIul(:,i) = ulnorm(i)*PHI(:,i);
end
mul = diag(PHIul'*M*PHIul);
kul = diag(PHIul'*K*PHIul);
figure; 
subplot(2,1,1); hold on;
plot(mul,'rx-');
xlabel('Mode #');
ylabel('Generalized Mass');
xlim([1 22]);
grid on;
title('Generalized Mass vs Mode #');
subplot(2,1,2); hold on;
plot(kul,'bx-');
xlabel('Mode #');
ylabel('Generalized Stiffness');
xlim([1 22]);
grid on;
title('Generalized Stiffness vs Mode #');
%%
%1d Mass Normalize. PHI'*M*PHI=I. Plot first two elastic modes.
%The shapes between 1a and 1d are the same. However, 1d's mode shapes are
%scaled in such a way as to make the modal mass matrix be identity.
mnorm = zeros(22,1);
PHImn = zeros(22,22);
for i = 1:22
   mnorm(i) = 1/(sqrt(PHI(:,i)'*M*PHI(:,i)));
   PHImn(:,i) = mnorm(i)*PHI(:,i);
end

figure; 
subplot(2,1,1);hold on;
plot(PHImn(:,3),'r-'); 
ylabel('Amplitude');
xlabel('Mode #');
xlim([1 22]);
grid on;
title('Mode Shape - 1st Elastic Mode');
subplot(2,1,2); hold on;
plot(PHImn(:,4),'b-');
ylabel('Amplitude');
xlabel('Mode #');
xlim([1 22]);
grid on;
title('Mode Shape - 2nd Elastic Mode');



