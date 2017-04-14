function [K,M,DOFind,DOFCOMP,Ksort,Msort] = getModalTAM(K,M,DOF,DOFred,...
PHI,PHIT)
% https://ay16-17.moodle.wisc.edu/prod/pluginfile.php/171625/mod_resource
%/content/2/Topic%209%20-%20Modal%20Reduction%20Pres.pdf

%Performs a Static Reduction of a model. Returns the reduced K and M
%matrices, the matrices are sorted in ascending node order.
%Input: K, M - Original mass and stiffness
% DOF - List of DOF (eg, 101.1, 101.2, 101.3 (node 101, x,y and z). This
% could also be a list of DOF indexes. Both DOF and DOFred need to have the
% same format.
% DOFred - List of DOF that should be kept.
% PHI - Provide the modes used for reduction.
% PHIT = index of target DOF.

if nargin ~= 6 %If PHI not provided, compute it.
    [PHI,~,~,~,sortindx] = getEigSort(K,M); 
end

nDOFred = length(DOFred); %# kept dof.
nDOF = length(DOF); % # dof total.
DOFind = zeros(nDOFred,1);
for i = 1:nDOFred
    [tmpind,~] = find(DOF==DOFred(i));
    DOFind(i) = tmpind;
end
DOFCOMP = setxor([1:nDOF],DOFind);

Ksort = [K(DOFind,DOFind) K(DOFind,DOFCOMP);...
    K(DOFCOMP,DOFind) K(DOFCOMP, DOFCOMP)];
Msort = [M(DOFind,DOFind) M(DOFind,DOFCOMP);...
    M(DOFCOMP,DOFind) M(DOFCOMP, DOFCOMP)];

%MODAL TAM
PHITs = PHI(DOFCOMP,PHIT);
PHITm = PHI(DOFind,PHIT);
T = [eye(nDOFred); PHITs*inv(PHITm'*PHITm)*PHITm'];
K = T'*Ksort*T;
M = T'*Msort*T;

end