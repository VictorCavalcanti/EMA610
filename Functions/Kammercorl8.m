function [C,p1,p2]=Kammercorl8(phi1,M,phi2)
%
%  This function normalizes the mode shapes wrt the
%  weighting matrix M, and then computes the triple
%  matrix multiplication 
%
%             C = phi1' * M * phi2
%
%  Use:  [C]=corl8(phi1,M,phi2);
%
%  Normalize Modes
%  ---------------
g = phi1'*M*phi1;
g = diag(g);
g = g.^(-1/2);
g=diag(g);
phi1 = phi1*g;
%
g = phi2'*M*phi2;
g = diag(g);
g = g.^(-1/2);
g=diag(g);
phi2 = phi2*g;
%
p1=phi1;
p2=phi2;
%
% Orthogonality Computation
% -------------------------
figure;
C = phi1'*M*phi2;
imagesc(abs(C));
colormap(parula(32)); colorbar;
axis equal;
xlim([0.5, numel(phi2(1,:))+0.5]);
ylim([0.5, numel(phi1(:,1))+0.5]);
set(gca,'XTick',1:numel(phi2(1,:)));
set(gca,'YTick',1:numel(phi1(:,1)));
%
end