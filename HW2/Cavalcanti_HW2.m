clear all; clc; close all;
load('station.mat');
%PICKED DOF FOR REDUCTION
DOF25 = [101.2 110.2 122.2 203.1  117.2 206.1 213.1 216.1 216.2 303.1 ...
    306.1 306.2 313.1 316.1 316.2 101.1 106.2 110.1 111.2 206.2 122.1 ...
    203.2 213.2 303.2 313.2]';

DOF15 = [101.2 110.2 122.2 203.1  117.2 206.1 213.1 216.1 216.2 303.1 ...
    306.1 306.2 313.1 316.1 316.2]';
DOF25 = sort(DOF25,'ascend');
DOF15 = sort(DOF15,'ascend');


%% HARD CODING, IMPLEMENTED FUNCTION FOR STATIC TAM INSTEAD
% DOF25ind = zeros(25,1);
% DOF15ind = zeros(15,1);
% for i = 1:25
%     [tmpind,~] = find(DOF==DOF25(i));
%     DOF25ind(i) = tmpind;
%     if i<=15
%        [tmpind,~] = find(DOF==DOF15(i));
%        DOF15ind(i) = tmpind;
%     end
% end
% DOF25COMP = setxor([1:92],DOF25ind);
% DOF15COMP = setxor([1:92],DOF15ind);
% 
% K25sort = [K(DOF25ind,DOF25ind) K(DOF25ind,DOF25COMP);...
%     K(DOF25COMP,DOF25ind) K(DOF25COMP, DOF25COMP)];
% M25sort = [M(DOF25ind,DOF25ind) M(DOF25ind,DOF25COMP);...
%     M(DOF25COMP,DOF25ind) M(DOF25COMP, DOF25COMP)];
% 
% K15sort = [K(DOF15ind,DOF15ind) K(DOF15ind,DOF15COMP);...
%     K(DOF15COMP,DOF15ind) K(DOF15COMP, DOF15COMP)];
% M15sort = [M(DOF15ind,DOF15ind) M(DOF15ind,DOF15COMP);...
%     M(DOF15COMP,DOF15ind) M(DOF15COMP, DOF15COMP)];
% 
% T25 = [eye(25); -inv(K(DOF25COMP,DOF25COMP))*K(DOF25COMP,DOF25ind)];
% T15 = [eye(15); -inv(K(DOF15COMP,DOF15COMP))*K(DOF15COMP,DOF15ind)];
% 
% K25 = T25'*K25sort*T25;
% M25 = T25'*M25sort*T25;
% 
% K15 = T15'*K15sort*T15;
% M15 = T15'*M15sort*T15;
%%


%STATIC REDUCTION VIA FUNCTION
[K25,M25,DOFind25,DOFcomp25,~,~] = getStaticTAM(K,M,DOF,DOF25);
[K15,M15,DOFind15,DOFcomp15,K15sort,M15sort] = getStaticTAM(K,M,DOF,DOF15);

%EIGENPROBLEMS
[phi25,lam25] = eig(K25,M25);
wn25 = sqrt(abs(diag(lam25)));
[wn25,wn25ind] = sort(wn25,'ascend');
phi25 = phi25(:,wn25ind);
wn25hz = wn25./(2*pi);

[phi15, lam15] = eig(K15,M15);
wn15 = sqrt(abs(diag(lam15)));
[wn15,wn15ind] = sort(wn15,'ascend');
phi15 = phi15(:,wn15ind);
wn15hz = wn15./(2*pi);

% %MASS NORMALIZE MODES
% for i = 1:25
%  mnorm(i) = 1/(sqrt(phi25(:,i)'*M25*phi25(:,i)));
%  phi25(:,i) = mnorm(i)*phi25(:,i);
% end
% for i = 1:15
%  mnorm(i) = 1/(sqrt(phi15(:,i)'*M15*phi15(:,i)));
%  phi15(:,i) = mnorm(i)*phi15(:,i);
% end

%%
%TAM-FEM COMPARISON
%Get the mode shapes (exclude rigid body) 
%for all the DOF in the 25 and 15 DOF TAM
phi25 = phi25(:,4:end); %Remove rigid body modes
phi15 = phi15(:,4:end);
PHIfem25 = PHI(DOFind25,4:end); 
PHIfem15 = PHI(DOFind15,4:end);
 M25corl8 = M(DOFind25,DOFind25);
 M15corl8 = M(DOFind15,DOFind15);
%MASS NORMALIZE MODES TO ORIGINAL M
for i = 1:22
 mnorm(i) = 1/(sqrt(phi25(:,i)'*M25corl8*phi25(:,i)));
 phi25(:,i) = mnorm(i)*phi25(:,i);
end
for i = 1:12
 mnorm(i) = 1/(sqrt(phi15(:,i)'*M15corl8*phi15(:,i)));
 phi15(:,i) = mnorm(i)*phi15(:,i);
end

% SO25 = corl8(phi25,phi25,M25);
% CO25 = corl8(phi25, PHIfem25,M25);
% SO15 = corl8(phi15,phi15,M15);
% CO15 = corl8(phi15,PHIfem15,M15);

SO25 = corl8(phi25,phi25,M25corl8);
CO25 = corl8(phi25, PHIfem25,M25corl8);
SO15 = corl8(phi15,phi15,M15corl8);
CO15 = corl8(phi15,PHIfem15,M15corl8);

label25 = 4:2:24;
label15 = 4:2:14;

%Cross orthogonality for 25DOF static reduction
figure;
imagesc(abs(CO25));
colormap(parula(32)); colorbar;
caxis([0 1]);
axis equal;
ylim([0.5, 22.5]);
xlim([0.5,27.5]);
set(gca,'YTick',1:2:22);
set(gca, 'YTickLabel', label25); % Change y-axis ticks labels.
set(gca,'XTick',0:5:30);
xlabel('FEM Mode #');
ylabel('TAM Mode #'); 
title('Cross-Orthogonality 25 dof TAM-FEM');


% labels = [100 200 400 1000 2000 5000 10000 20000 50000];
% plot(y);
% set(gca, 'YTick', 1:length(labels)); % Change x-axis ticks



%Self orthogonoality for 25DOF static reduction
figure;
imagesc(abs(SO25));
colormap(parula(32)); colorbar;
caxis([0 1]);
axis equal;
ylim([0.5, 22.5]);
xlim([0.5,22.5]);
set(gca,'YTick',1:2:22);
set(gca, 'YTickLabel', label25); % Change y-axis ticks labels.
set(gca,'XTick',0:5:25);
xlabel('FEM Mode #');
ylabel('TAM Mode #');
title('Self-Orthogonality 25 dof TAM');
%Cross-O for 15DOF
figure;
imagesc(abs(CO15));
colormap(parula(32)); colorbar;
caxis([0 1]);
axis equal;
ylim([0.5, 12.5]);
xlim([0.5,27.5]);
set(gca,'YTick',1:2:12);
set(gca, 'YTickLabel', label15); % Change y-axis ticks labels.

set(gca,'XTick',0:5:30);
xlabel('FEM Mode #');
ylabel('TAM Mode #');
title('Cross-Orthogonality 15 dof TAM-FEM');
%Self orthogonoality for 15DOF static reduction
figure;
imagesc(abs(SO15));
colormap(parula(32)); colorbar;
caxis([0 1]);
axis equal;
ylim([0.5, 12.5]);
xlim([0.5,12.5]);
set(gca,'YTick',1:2:12);
set(gca, 'YTickLabel', label15); % Change y-axis ticks labels.

set(gca,'XTick',0:5:15);
xlabel('FEM Mode #');
ylabel('TAM Mode #');
title('Self-Orthogonality 15 dof TAM-FEM');




%%
%LSIM WITH FORCE VECTOR PROVIDED
% https://ay16-17.moodle.wisc.edu/prod/pluginfile.php/118403/ ...
%  mod_resource/content/1/Matlab%20Simulation.pdf
%15DOF TAM- First 11 modes
%Setting up force vectors, Force is applied at node 117.2.
F15 = [zeros(2,1000); F'; zeros(12,1000)];

phi11 = phi15(:,1:11);
m = phi11'*M15*phi11; %get modal mass and stiffness with 11 modes
k = phi11'*K15*phi11;
A11 = [zeros(11), eye(11); -m\k zeros(11)]; %Set up state space
B11 = [zeros(11);inv(m)];
G11 = [ -m\k zeros(11)];
D11 = inv(m);
sys11 = ss(A11,B11,G11,D11); %Set up sys
% u = (m\(phi11'))*F15;
u = phi11'*F15; %Set up input channel vector
x0 = zeros(22,1);
t = 0:0.01:(10-0.01);
[y] = lsim(sys11,u',t,x0);
yphys = phi11*y';

% FEM LSIM
kfem = PHI'*K*PHI;
mfem = PHI'*M*PHI;
ufem = PHI'*[zeros(33,1000); F'; zeros(58,1000)];
A = [zeros(30), eye(30); -mfem\kfem zeros(30)];
B = [zeros(30); inv(mfem)];
G = [-mfem\kfem zeros(30)];
D = inv(mfem);
sysfem = ss(A,B,G,D);
x0fem = zeros(60,1);
[yfem] = lsim(sysfem,ufem',t,x0fem);
yfemphys = PHI*yfem';


%MODAL TAM
% phitam = PHI([DOFind15, DOFcomp15],:);
PHITs = PHI(DOFcomp15,1:11);
PHITm = PHI(DOFind15,1:11);
TTam = [eye(15); PHITs*inv(PHITm'*PHITm)*PHITm'];
KT = TTam'*K15sort*TTam;
MT = TTam'*M15sort*TTam;
[phiT,lamT] = eig(KT,MT);
wnT = sqrt(abs(diag(lamT)));
[wnT,wnTind] = sort(wnT,'ascend');
phiT = phiT(:,wnTind);
wnThz = wnT./(2*pi);


%LSIM WITH MODAL TAM
kT = phiT(:,1:11)'*KT*phiT(:,1:11); %Only use first 11 modes from modal TAM
mT = phiT(:,1:11)'*MT*phiT(:,1:11);
uT = phiT(:,1:11)'*[zeros(2,1000); F'; zeros(12,1000)];
At = [zeros(11), eye(11); -mT\kT zeros(11)];
Bt = [zeros(11); inv(mT)];
Gt = [-mT\kT zeros(11)];
Dt = inv(mT);
sysT = ss(At,Bt,Gt,Dt);
x0t = zeros(22,1);
[yT] = lsim(sysT,uT',t,x0t);
yTphys = phiT(:,1:11)*yT';


%PLOTTING ALL LSIM RESULTS
figure;
plot(t,yfemphys(55,:),'k--','LineWidth',1.5); hold on;
plot(t, yphys(6,:),'b','LineWidth',1);
plot(t,yTphys(6,:),'r','LineWidth',1);
xlabel('Time (s)');
ylabel('Physical Displacement (units not provided)');
legend('FEM', '15DOF Static TAM', '11 Mode 15DOF Modal TAM');
title('Physical Acceleration of DOF 206x due to loading at 117y');

%PLOTTING FREQUENCY ERRORS
%Novel way doesn't seem very good here. Tabule results instead.
% figure; hold on;
% plot(w(1:11),wn15hz(1:11),'k-x');
% plot(w(1:11),wn25hz(1:11),'r-o');
% plot(w(1:11),wnThz(1:11),'b-o');