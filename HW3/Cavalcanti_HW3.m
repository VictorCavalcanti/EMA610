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
I = find(E>=0.047); 
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

[Efi,EFIdof,KEdof] = getEffectiveIndependence(PHIK,MF,KF,numel(modindx)+5,dofO);


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
disp('The DOF set selected by EFI is:');
disp(DOF1(EFIdof)');
disp('The DOF set selected by KE is:');
disp(DOF1(KEdof)');


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

% [efiPHI,efiwn] = eig(Kefi,Mefi);
% efiwnhrz = sqrt(abs(diag(efiwn)));
% [efiwnhrz, efiascindx] = sort(efiwnhrz,'ascend'); %Sort frequencies, ascend.
% efiPHI = efiPHI(:,efiascindx); %sort shapes based on freqs.
% 
% [kePHI,kewn] = eig(Kke,Mke);
% kewnhrz = sqrt(abs(diag(kewn)));
% [kewnhrz, keascindx] = sort(kewnhrz,'ascend'); %Sort frequencies, ascend.
% kePHI = kePHI(:,keascindx); %sort shapes based on freqs.
PHImacefi = PHI(efiREDind,:);
PHImacke = PHI(keREDind,:);
%MAC
efimac = mac(efiPHI,PHImacefi);
title('Modal Assurance Criterion (MAC) between PHI_e_f_i and PHI');
ylabel('PHI_e_f_i Mode #');
xlabel('PHI Mode #');
kemac = mac(kePHI,PHImacke);
title('Modal Assurance Criterion (MAC) between PHI_k_e and PHI');
ylabel('PHI_k_e Mode #');
xlabel('PHI Mode #');

%Cross and Self Orthogonalization

%Frequency error
