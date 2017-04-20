close all; clear all; clc; 
load('gpsc.mat');

wtmass = 386.4; %parameter to convert weight to mass for M.
% M = wtmass.*M; %Convert to propper units.

%Sort DOF1 and DOF2 in ascending order. Reorder M and K according to DOF1
%sort, sort PHI according to DOF2 sort. 
[DOF1, DOF1ind] = sort(DOF1,'ascend');
[DOF2, DOF2ind] = sort(DOF2,'ascend');
M = M(DOF1ind,DOF1ind);
K = K(DOF1ind,DOF1ind);
PHI = PHI(DOF2ind,:);
PHIFI = PHI; %fixed interface modes;

%interface dof set
aSET = [44.1; 44.2; 44.3; 45.1; 45.2; 45.3; 48.1; 48.2; 48.3; 49.1; ...
    49.2; 49.3];

%Compute RBM about node 50
PHIR = [  -inv(K(1:end-6,1:end-6))*K(1:end-6,end-5:end); eye(6)];
PHIR = PHIR(1:end-3,:); %Remove the rotational DOF from PHIR
DOF1 = DOF1(1:end-3); %Remove the rotation DOF from list.



%Static reduction of the 3 rotations at DOF 50 (50.4,50.5,50.6);
T = [eye(150); -inv(K(151:end,151:end))*K(151:end,1:150)];
M = T'*M*T;
K = T'*K*T;



%% Generate a function that computes effective mass and use it to rank the
%fixed interface modes of the craft (PHIsort). Select a set of target modes
%based on the fact that each of them needs at least 5% of the effective
%mass in each of the six rigid body directions. Determine the total
%effective mass for the target mode set in each rigid body direction and
%comment on the sets dynamic completeness.



% M = wtmass.*M; %Convert to propper units.
[E,MR,MRo,Msort,Ksort,dofA,dofO] = getEffectiveMass(M,K,PHIFI,aSET,DOF1,PHIR);
Esum = sum(E);
%%
% DOF1 = dofO;
% MR = wtmass.*MR;
% MRo = wtmass.*MRo;
%%
disp('Effective Mass (30 modes):');
disp(Esum);

%Selects all modes that contribute to more than 4.65% of mass in any rigid
%body direction. The assignment specifies 5%, however Dr. Kammer suggested
%"Things close to 5% RBM".
I = find(E>=0.048); 
modindx = sort(mod(I,numel(E(:,1))),'ascend');
modindx = unique(modindx);
masscmplt = sum(E(modindx,:));
disp('Selected Modes Total Effective Mass:');
disp(masscmplt);
disp('The Modes that should be kept are:');
disp(modindx');

%Fixed interface M and K. These matrices were sorted in the
%getEffectiveMass call. (O-set partition)
PHIK = PHIFI(:,modindx); %Kept fixed-interface modes. 
MF = Msort(13:end,13:end); 
KF = Ksort(13:end,13:end);  



%Generate functions to rank candidate sensor locations based on Modal
%Kinetic Energy and Effective Independence. Pick a single initial candidate
%set of sensor locations and use these functions to select a final set that
%will identify your target modes. The final sensor location should have
%n-target modes + 5 sensors.

[Efi,EFIdof,KEdof,efiDetQ,keDetQ] = getEffectiveIndependence(PHIK,MF,KF,numel(modindx)+5,dofO);


%Use MAC and any other measure to determine which of the sensor sets
%(initial and final) produces the most independent target mode partitions
%and the greatest target mode response signal strength.

EFIdof = sort(EFIdof,'ascend'); %DOF in original 150x150 mass matrix.
KEdof = sort(KEdof,'ascend');

for i = 1:numel(EFIdof)                        
    [efiDOFind(i),~] = find(EFIdof(i)==dofO);
end
for i = 1:numel(KEdof)                        
    [keDOFind(i),~] = find(KEdof(i)==dofO);
end
efiDOFind = efiDOFind';
keDOFind = keDOFind';
disp('The DOF set selected by EFI are:');
disp(DOF1(EFIdof)');
disp('The DOF set selected by KE are:');
disp(DOF1(KEdof)');
grid on;

%% DEBUG ONLY
% i = [ 130:135, 142:147]';
% j = setxor(1:150,i);
% [PHI,D,wn,wnhz,sortindx] = getEigSort(K(j,j),M(j,j));
%%


%Static reductions to EFIdof and KEdof.
[Kke,Mke,keREDind,keREDcomp,Kkesort,Mkesort] = ...
getStaticTAM(KF,MF,[1:138]',keDOFind);
[Kefi,Mefi,efiREDind,efiREDcomp,Kefisort,Mefisort] = ...
getStaticTAM(KF,MF,[1:138]',efiDOFind);
%Eigenproblem for EFIdof and KEdof.
[efiPHI,efiD,efiwn,efiwnhz,efisortindx] = getEigSort(Kefi,Mefi);
[kePHI,keD,kewn,kewnhz,kesortindx] = getEigSort(Kke,Mke);
PHImacefi = PHI(efiREDind,modindx);
PHImacke = PHI(keREDind,modindx);

%MAC
efimac = mac(efiPHI,PHImacefi);
set(gca, 'XTickLabel', modindx);
title('Modal Assurance Criterion (MAC) between PHI_e_f_i and PHI');
xlabel('Target Fixed-Interface Mode #');

efiSmac = mac(efiPHI,efiPHI);
title('Self - Modal Assurance Criterion (MAC) for PHI_e_f_i');
ylabel('PHI_e_f_i Mode #');
xlabel('PHI_e_f_i Mode #');

kemac = mac(kePHI,PHImacke);
set(gca, 'XTickLabel', modindx);
title('Modal Assurance Criterion (MAC) between PHI_k_e and PHI');
xlabel('Target Fixed-Interface Mode #');

keSmac = mac(kePHI,kePHI);
title('Self - Modal Assurance Criterion (MAC) for PHI_k_e');
ylabel('PHI_k_e Mode #');
xlabel('PHI_k_e Mode #');

%Cross and Self Orthogonality
% [keCO] = corl8(kePHI,PHImacke,Mke);
[Cke,~,~] = Kammercorl8(kePHI,Mke,PHImacke);
set(gca, 'XTickLabel', modindx);
title('Cross Orthogonality (CO) between PHI_k_e and PHI');
ylabel('PHI_k_e Mode #');
xlabel('Target Fixed-Interface Mode #');

% [keSO] = corl8(kePHI,kePHI,Mke);
[Sefi,~,~] = Kammercorl8(kePHI,Mke,kePHI);
% set(gca, 'XTickLabel', modindx);
title('Self Orthogonality (SO) for PHI_k_e');
ylabel('PHI_k_e Mode #');
xlabel('Target Fixed-Interface Mode #');

% [efiCO] = corl8(efiPHI,PHImacefi,Mefi);
[Cefi,p1,p2] = Kammercorl8(efiPHI,Mefi, PHImacefi);
set(gca, 'XTickLabel', modindx);

title('Cross Orthogonality (CO) between PHI_e_f_i and PHI');
ylabel('PHI_e_f_i Mode #');
xlabel('Target Fixed-Interface Mode #');

% [efiSO] = corl8(efiPHI,efiPHI,Mefi);
[Sefi,~,~] = Kammercorl8(efiPHI,Mefi, efiPHI);
% set(gca, 'XTickLabel', modindx);
title('Self Orthogonality (SO) for PHI_e_f_i');
ylabel('PHI_e_f_i Mode #');
xlabel('Target Fixed-Interface Mode #');

%Frequency error
% efimodematching = [1,1; 2,4; 3,3; 4,5; 5,6; 7,8; 10,9; 11,12; 12,13; 15,18; 16,19]; %FULL SET
kemodematching = [1,1; 2,3; 3,4; 4,5; 5,6; 7,8; 8,9; 11,12; 12,13; 14,18; 15,19; 16,20; 17,21; 18,22]; %FULL SET (REDUCED IS THE SAME)

efimodematching = [1,1; 2,3; 3,4; 4,6; 5,5; 6,8; 8,9; 10,12; 11,13; 13,18; 14,19; 15,21]; %REDUCED BEGINING SET

% kemodematching = [1,1; 2,3; 3,4; 4,5; 5,6; 7,9; 11,12; 12,13; 14,18; 15,19;16,20;17,21; 18,22];
eficompare = [w(efimodematching(:,2)) efiwnhz(efimodematching(:,1))];
kecompare = [w(kemodematching(:,2)) kewnhz(kemodematching(:,1))];