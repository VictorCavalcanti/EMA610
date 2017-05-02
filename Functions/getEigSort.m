function [PHI,D,wn,wnfreq,sortindx] = getEigSort(K,M)
%Functions that returns sorted eigenvalues, eigenvectors (modes),
%frequencies in radians and hertz
if nargin == 1
    [PHI,D] = eig(K);
else
    [PHI,D] = eig(K,M);
end
D = diag(D);
wn = sqrt(abs(D));
[wn, sortindx] = sort(wn,'ascend'); %Sort frequencies, ascend.
PHI = PHI(:,sortindx); %sort shapes based on freqs.
%Mass Normalize Modes
g = PHI'*M*PHI;
g = diag(g);
g = g.^(-1/2);
g=diag(g);
PHI = PHI*g;
%

D = D(sortindx);
wnfreq = wn./(2*pi);

end