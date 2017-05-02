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
title('Cross Orthogonality (CO) between PHI_s_t_a_t_i_c and PHI_F_E_M');
ylabel('PHI_s_t_a_t_i_c Mode #');
xlabel('Target Mode #');

%%
%MODAL TAM
[Km,Mm,DOFindm,DOFCOMPm,Ksortm,Msortm,Tmod] = getModalTAM(K,M,DOFindx,ASETindx,...
    PHI,TMODindx);
%Solve eigenvalue problem and sort modes ascending.
[PHIm,mD,mwn,mwnhz,msortindx] = getEigSort(Km,Mm);
%Plot bits
[COm,~,~] = Kammercorl8(PHIm,M(ASETindx,ASETindx), PHIT(ASETindx,:));
set(gca, 'XTickLabel', TMODindx);
title('Cross Orthogonality (CO) between PHI_m_o_d_a_l and PHI_F_E_M');
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

%Extract Mode Shapes
mindx = [ 19, 32, 48, 68, 92];

figure;
% subplot(2,1,1); 
hold on;
subplot(3,1,1); hold on;
plot(w,(360/pi)*angle(g(1,:)),'LineWidth',1.5);
% plot(w,(360/pi)*angle(g(2,:)),'r','LineWidth',1.5);
ylim([-360 360]);
xlim([0 512]);
ystart = -360*ones(1,5);
yend = 360*ones(1,5);
for idx = 1 : numel(ystart)
    plot([2.*mindx(idx) 2.*mindx(idx)], [ystart(idx) yend(idx)],...
        'k--','Linewidth',1.25);
end
% plot(w,(360/pi)*angle(g(3,:)),'k');
% plot(w,(360/pi)*angle(g(4,:)),'c');
% plot(w,(360/pi)*angle(g(5,:)),'y');

title('Phase Angle of drive point FRF');
ylabel('Angle (Degrees)');
xlabel('Frequency (Hz)');
% grid on;

% figure; 
subplot(3,1,2);

% plot(w,real(g(1,:))); 
% Inertance: A/F
semilogy(w,abs(g(1,:)),'LineWidth',1.5); hold on;
% semilogy(w,abs(g(2,:)),'r','LineWidth',1.5);
% axis([0 512 -s s ]);
xlim([0 512]);
% ystart = -s.*ones(1,5);
ystart = 0.01.*ones(1,5);
yend = 100.*ones(1,5);
for idx = 1 : numel(ystart)
    plot([2.*mindx(idx) 2.*mindx(idx)], [ystart(idx) yend(idx)],...
        'k--','Linewidth',1.25);
end
% semilogy(w,abs(g(3,:)),'k');
% semilogy(w,abs(g(4,:)),'c');
% semilogy(w,abs(g(4,:)),'c');
title('Inertance vs Frequency of drive point FRF');
ylabel('Magnitude of Inertance');
xlabel('Frequency (Hz)');

% grid on;

%Compliance: X/F
h = g./repmat(w,1,5)';
% figure; 
subplot(3,1,3); 

semilogy(w,abs(h(1,:)),'LineWidth',1.5); hold on;
% semilogy(w,abs(h(2,:)),'r','LineWidth',1.5);
xlim([0 512]);

ystart = 0.0001.*ones(1,5);
yend = 10.*ones(1,5);
for idx = 1 : numel(ystart)
    plot([2.*mindx(idx) 2.*mindx(idx)], [ystart(idx) yend(idx)],...
        'k--','Linewidth',1.25);
end
% semilogy(w,abs(h(3,:)),'k');
% semilogy(w,abs(h(4,:)),'c');
% grid on;
title('Compliance vs Frequency of drive point FRF');
ylabel('Magnitude of Compliance');
xlabel('Frequency (Hz)');


phic = g(:,mindx);
%Realize modes according to:
%   https://sem.org/wp-content/uploads/2016/07/
%   sem.org-IMAC-XI-11th-Int-11-39-6-Realization-Complex-Mode-Shapes.pdf
% Equation 2.
PHIt = real(phic)+imag(phic)*inv(real(phic).'*real(phic))*real(phic).'*imag(phic);
%Expand the test modes to full DOF size. (using modal transformation used in part 1)
PHIrFull = Tmod*PHIt;

%Cross-Orthogonality
[TAMTESTCO,~,~] = Kammercorl8(PHIm,M(ASETindx,ASETindx), PHIt);
%Plot bits
set(gca, 'XTickLabel', TMODindx);
title('Cross Orthogonality (CO) between PHI_T_A_M and PHI_T_E_S_T');
ylabel('TAM Mode #');
xlabel('Target Mode #');
%MAC
mac = mac(PHIm,PHIt);
title('Cross - MAC for PHI_T_A_M vs PHI_T_E_S_T');
ylabel('TAM Mode #');
xlabel('TEST Mode #');

%COMPARE FEM/TEST TARGET MODE CORRELATION
[FEMTESTCO,~,~] = Kammercorl8(PHI([ASETindx ;DOFCOMPm],TMODindx),Msortm,PHIrFull );
%Plot bits
set(gca, 'XTickLabel', TMODindx);
set(gca, 'YTickLabel', TMODindx);
title('Cross Orthogonality (CO) between PHI_F_E_M and PHI_T_E_S_T');
ylabel('FEM Mode #');
xlabel('TEST Mode #');
ylim([0.5, 5+0.5]);

[~,plotidx] = sort([ASETindx; DOFCOMPm],'ascend');

% figure; 
% subplot(5,1,1); 
% plot(PHI(:,4),'r','LineWidth',1.5); hold on;
% plot(5.*PHIrFull(plotidx,1),'b','LineWidth',1.5);
% 
% subplot(5,1,2); 
% plot(PHI(:,5),'r','LineWidth',1.5); hold on;
% plot(-5.*PHIrFull(plotidx,2),'b','LineWidth',1.5);
% 
% subplot(5,1,3); 
% plot(PHI(:,6),'r','LineWidth',1.5); hold on;
% plot(-2.*PHIrFull(plotidx,3),'b','LineWidth',1.5);
% 
% subplot(5,1,4); 
% plot(PHI(:,7),'r','LineWidth',1.5); hold on;
% plot(2.*PHIrFull(plotidx,4),'b','LineWidth',1.5);
% 
% subplot(5,1,5); 
% plot(PHI(:,8),'r','LineWidth',1.5); hold on;
% plot(-5.*PHIrFull(plotidx,5),'b','LineWidth',1.5);

figure; 
subplot(4,1,1); 
plot(PHI(:,8),'r','LineWidth',1.5); hold on;
plot(-2.*PHIrFull(plotidx,4),'b','LineWidth',1.5);
title('FEM mode 5 & Test Mode 4');
subplot(4,1,2); 
plot(PHI(:,6),'r','LineWidth',1.5); hold on;
plot(-5.*PHIrFull(plotidx,2),'b','LineWidth',1.5);
title('FEM mode 6 & Test Mode 5');
subplot(4,1,3); 
plot(PHI(:,7),'r','LineWidth',1.5); hold on;
plot(2.*PHIrFull(plotidx,3),'b','LineWidth',1.5);
title('FEM mode 7 & Test Mode 6');

subplot(4,1,4); 
plot(PHI(:,8),'r','LineWidth',1.5); hold on;
plot(-2.*PHIrFull(plotidx,4),'b','LineWidth',1.5);
title('FEM mode 8 & Test Mode 7');

legend('FEM Mode', 'Test Mode');