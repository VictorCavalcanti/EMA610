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
% nt = size(a1,1)/2;
% [Y1,w1] = fft_easy(a1(:,1),ts);
% [Y2,~] = fft_easy(a2(:,1),ts);
% [Y3,~] = fft_easy(a3(:,1),ts);
% Y = [Y1, Y2, Y3];
y = zeros(size(a1,1),size(a1,2),ntest);
f = zeros(size(f1,1),size(f1,2),ntest);
% y = Y;
% f = F;

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
% CS = zeros(ns,na);
% AS = zeros(na,na);
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
%     pxy = Y*F';
%     pxx = F*F';
%     cs = CS+Y(:,i)*U(:,i)';
%     as = AS+U(:,i)*U(:,i)';
Pxy=Pxy+pxy;
Pxx=Pxx+pxx;
end
% g = Pxy*inv(Pxx);

for j=1:nt                         % loop over data points
  g(:,(j-1)*na+ia)=Pxy(:,(j-1)*na+ia)*...
  inv(Pxx(:,(j-1)*na+ia));         % build frf matrix
end
figure;
subplot(2,1,1); 
plot(w,angle(g(1,:)));
subplot(2,1,2);
plot(w,real(g(1,:)));

% g = AS\CS;

%  Build frequency response function matrix
%  ----------------------------------------
% ia=(1:na)';                        % input counter
% is=(1:ns)';                        % output counter
% Pxy=zeros(ns,nt*na);
% Pxx=zeros(na,nt*na);
% %
% for i = 1:nave                     % loop over tests
%     uf=fft(u);                       % DFT of input
%     yf=fft(y);                       % DFT of output
%     uf=uf.';
%     yf=yf.';
%     %
%     for j=1:nt                       % loop over data points
%         xy=yf(:,j)*uf(:,j)';           % cross spectral density at j
%         xx=uf(:,j)*uf(:,j)';           % spectral desity at j
%         pxy(:,(j-1)*na+ia)=xy;         % build csd matrix
%         pxx(:,(j-1)*na+ia)=xx;         % build sd matrix
%     end
%     %
%     Pxy=Pxy+pxy;                     % cross spectral density
%     Pxx=Pxx+pxx;                     % spectral density
% end
%