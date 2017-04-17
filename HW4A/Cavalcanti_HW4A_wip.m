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
ntest = 3;
nt = 257;

y = zeros(size(a1,1),size(a1,2),ntest);
f = zeros(size(f1,1),size(f1,2),ntest);

y(:,:,1) = a1;
y(:,:,2) = a2;
y(:,:,3) = a3;
f(:,:,1) = f1;
f(:,:,2) = f2;
f(:,:,3) = f3;

ns = size(y,2);
na = size(f,2);
ia=(1:na)';                        % input counter
is=(1:ns)';                        % output counter
Pxy=zeros(ns,nt*na);
Pxx=zeros(na,nt*na);

for i = 1:ntest %loop over tests
    [Y,w] = fft_easy(y(:,:,i),ts);
    [F,~] = fft_easy(f(:,:,i),ts);
    Y = Y.';
    F = F.';
    for j = 1:nt
        xy = Y(:,j)*F(:,j)';
        xx = F(:,j)*F(:,j)';
        pxy(:,(j-1)*na+ia)=xy;
        pxx(:,(j-1)*na+ia)=xx;
    end
Pxy=Pxy+pxy;
Pxx=Pxx+pxx;
end

for j=1:nt                         % loop over data points
  g(:,(j-1)*na+ia)=Pxy(:,(j-1)*na+ia)*...
  inv(Pxx(:,(j-1)*na+ia));         % build frf matrix
end
w = w./(2*pi); %fft easy returns frequencies in radians.
figure;
subplot(2,1,1); 
plot(w,angle(g(1,:))); hold on;
plot(w,imag(g(1,:)),'r');
subplot(2,1,2);
plot(w,real(g(1,:))); hold on;
plot(w,abs(g(1,:)),'r');
