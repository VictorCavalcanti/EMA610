function [Efi,DOF,xtraDOF] = getEffectiveIndependence(PHI,M,K, ntargetdof,DOF)
%Effective Independce can be explained here:
%  https://ay16-17.moodle.wisc.edu/prod/pluginfile.php/171634/
%  mod_resource/content/3/eCOWI_Resources/lecturepre/
%  Topic%2010.2%20-%20Sensor%20Placement%20using%20Effective
%  %20Independence%20Pres.pdf

%And here:
%  https://ay16-17.moodle.wisc.edu/prod/pluginfile.php/171635/
%  mod_resource/content/1/AIAA-JGCD-91.pdf


%It uses Modal Kinetic Energy to pick an initial candidate set of DOF which
%is 2x ntargetdof. It then iteratively removes dof by ranking via effective
%independence until a set of target dof is chosen (with ntargetdof dof).



%xtraDOF - DOF FOR HW3. SET OF nM+5 dof based solely on MKE.

[n,nM] = size(PHI);

for i = 1:n
    for j = 1:nM
        KE(i,j) = PHI(i,j)*M(i,:)*PHI(:,j);
    end
end

KEnorm = (1/nM)*sum(KE,2);
%idx is the index of the positions in the PHI vector of the sorted kinetic
%energy by modes.
[KEsort, idx] = sort(KEnorm,'descend');
xtraDOF = DOF(idx(1:19));
%DOF give the indexes of the original DOF list
DOF = DOF(idx(1:2*ntargetdof));

%Reduce our target modes to the initial set of DOF selected by Modal
%Kinetic Energy
PHI = PHI(idx(1:2*ntargetdof),:);
k = 2*ntargetdof;
while k > ntargetdof
    
    A = zeros(nM,nM);
    for i = 1:size(PHI,1)
        A = A+PHI(i,:)'*PHI(i,:);
    end
    [V,D] = eig(A);
    [D, ascindx] = sort(diag(D),'ascend'); %Sort eigenvalues.
    V = V(:,ascindx); %sort eigenvectors accordingly.
    G = (PHI*V).^2;
    % F = G.*repmat((1./D'),38,1);
    F = G*inv(diag(D));
    Efi = F*ones(nM,1);
    [Efi,index] = sort(Efi,'descend');
%Remove the lowest contributing dof from original DOF list and from PHI.
    DOF = DOF(index(1:end-1)); 
    PHI = PHI(index(1:end-1),:);
    
    k = k-1;
end
   
end

