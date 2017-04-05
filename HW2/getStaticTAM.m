function [K,M,DOFind,DOFCOMP,Ksort,Msort] = getStaticTAM(K,M,DOF,DOFred)
%https://ay16-17.moodle.wisc.edu/prod/pluginfile.php/171622/mod_resource...
%/content/3/Topic%208%20-%20Static%20Reduction%20Pres.pdf

%Performs a Static Reduction of a model. Returns the reduced K and M
%matrices, the matrices are sorted in ascending node order.
%Input: K, M - Original mass and stiffness
% DOF - List of DOF (eg, 101.1, 101.2, 101.3 (node 101, x,y and z)
% DOFred - List of DOF that should be kept (same format as above)

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
T = [eye(nDOFred); -inv(K(DOFCOMP,DOFCOMP))*K(DOFCOMP,DOFind)];
K = T'*Ksort*T;
M = T'*Msort*T;

end