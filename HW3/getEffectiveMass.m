function [E,MR, MRo,M,K] = getEffectiveMass(M,K,PHI, DOFint, DOF, PHIR )
% Effective Mass can be described here:
%  https://ay16-17.moodle.wisc.edu/prod/pluginfile.php/171644/...
%  mod_resource/content/1/Topic%2012.1%20-%20Selection%20of%20Target...
%  %20Modes%20-%20Effective%20Mass%20Pres.pdf

% Returns E, the effective mass matrix. It's columns sum to the diagonal
% terms of the rigid body mass matrix Mro (PHI'*Moo*PHI)
% Returns MR, the rigid body mass matrix, MRo, the o-set partition of MR.
% M and K, sorted mass and stiffness.

% Takens in the mass matrix M, FIXED-INTERFACE MODES PHIro,
% modes PHI, the list of DOF at the interface DOFint and
% the total list of DOF, DOF.
% Mass matrix and rows of PHI should be sorted according to DOF.

nDOFint = length(DOFint); % # kept dof.
nDOF = length(DOF); % # dof total.
DOFIND = zeros(nDOFint,1); % INDEX OF ASET wrt M and K
for i = 1:nDOFint
    [tmpind,~] = find(DOF==DOFint(i));
    DOFIND(i) = tmpind;
end
DOFCOMP = setxor(1:nDOF,DOFIND); % INDEX O SET DOF wrt M and K
DOF = DOF([DOFIND;DOFCOMP],:); % ordered DOF list, [ASET;OSET]
Moo = M(DOFCOMP,DOFCOMP);
for i = 1:numel(PHI(1,:))
 mnorm(i) = 1/(sqrt(PHI(:,i)'*Moo*PHI(:,i)));
 PHI(:,i) = mnorm(i)*PHI(:,i);
end
% Moa = M(DOFind,DOFCOMP);
% Mao = M(DOFCOMP,DOFind);
% Maa = M(DOFCOMP,DOFCOMP);
%M AND K PARTITIONED ACCORDING TO DOF 
M = [ M(DOFIND,DOFIND), M(DOFIND,DOFCOMP); ...
    M(DOFCOMP,DOFIND),M(DOFCOMP,DOFCOMP)];
K = [ K(DOFIND,DOFIND), K(DOFIND,DOFCOMP); ...
    K(DOFCOMP,DOFIND),K(DOFCOMP,DOFCOMP)];

%Orthogonality Check
surf(PHI'*Moo*PHI);
% OSET MASS WEIGHTED FIXED INTERFACE MODES ORTHOGONALITY CHECK.

%Compute PHIro *** NOTE, this restricts rigid body modes from the first 6
%DOF of the interface.

% PHIr = [ eye(6); -inv(K(7:end,7:end))*K(7:end,1:6)];
% PHIRo =  PHIR(nDOFint+1:end,:);
PHIRo = PHIR(DOFCOMP,:);
PHIR = PHIR([DOFIND; DOFCOMP],:);
MR = PHIR'*M*PHIR;
MRo = PHIRo'*Moo*PHIRo; %Modal mass matrix contributing to o-set partition.


%%
A = repmat((1./diag(MRo))',30,1);
E = ([PHI'*Moo*PHIRo].^2).*A;
% E = ([PHI'*Moo*PHIRo].^2)*(1./diag(MRo));
end