close all; clear all; clc;
%Victor Cavalcanti
%EMA610 - HW4A
%% 
%STATIC TAM
load('beam601.mat');
%Decipher what L is: hrz? radians?
l = diag(L);
lhz = abs(l.^(.5))./(2*pi);

TMODindx = 4:8; %Target mode indexes
ASETindx = [1;13;17;29;41]; %Target accelerometer locations
PHIT = PHI(:,TMODindx); %Target modes.
DOFindx = [1:numel(DOF)]'; %All dof indexes (original indexes)
%Get Static TAM
[Ks,Ms,DOFinds,DOFCOMPs,Ksorts,Msorts] = getStaticTAM(K,M,DOFindx,ASETindx); 


%Solve eigenvalue problem and sort modes ascending.
[PHIs,sD,swn,swnhz,ssortindx] = getEigSort(Ks,Ms); 
[COs,~,~] = Kammercorl8(PHIs,M(ASETindx,ASETindx), PHIT(ASETindx,:));
%Plot bits
set(gca, 'XTickLabel', TMODindx);
title('Cross Orthogonality (CO) between PHI_s_t_a_t_i_c and PHI_t_e_s_t'); 
ylabel('PHI_s_t_a_t_i_c Mode #');
xlabel('Target Mode #');

%%
%MODAL TAM
[Km,Mm,DOFindm,DOFCOMPm,Ksortm,Msortm] = getModalTAM(K,M,DOFindx,ASETindx,...
PHI,TMODindx);
%Solve eigenvalue problem and sort modes ascending.
[PHIm,mD,mwn,mwnhz,msortindx] = getEigSort(Km,Mm); 
%Plot bits
[COm,~,~] = Kammercorl8(PHIm,M(ASETindx,ASETindx), PHIT(ASETindx,:));
set(gca, 'XTickLabel', TMODindx);
title('Cross Orthogonality (CO) between PHI_m_o_d_a_l and PHI_t_e_s_t'); 
ylabel('PHI_m_o_d_a_l Mode #');
xlabel('Target Mode #');
%%
%GENERATE FREQUENCY RESPONSE FUNCTION FOR NODE 1 (DRIVE POINT), MAGNITUDE
%AND PHASE ANGLE.
%Sampling rate for taptest - 1024Hz.
%https://ay16-17.moodle.wisc.edu/prod/pluginfile.php/171656/
%mod_resource/content/5/eCOWI_Resources/lecturepre/
%Topic%2013%20-%20Basics%20of%20Modal%20Testing%20Pres.pdf

load('taptest.mat');

ts = 1/1024;
[Y1,w1] = fft_easy(a1,ts);
[Y2,w2] = fft_easy(a2,ts);
[Y3,w3] = fft_easy(a3,ts);

[U1,wu2] = fft_easy(f1,ts);
[U2,wu2] = fft_easy(f2,ts);
[U3,wu2] = fft_easy(f3,ts);

% fft_easy(a1(:,1),ts);
% phs = angle(fftshift(Xfft));
% figure;
% hold on;
% plot(phs,'r');
% plot(angle(Xfft),'b');