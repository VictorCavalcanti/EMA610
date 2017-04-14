close all; clear all; clc;
%Victor Cavalcanti
%EMA610 - HW4A
%% 
%STATIC TAM
load('beam601.mat');
%Decipher what L is: hrz? radians?
l = diag(L);
lhz = l./(2*pi);

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
load('taptest.mat');
