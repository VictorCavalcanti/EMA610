function [MAC] = mac(PHI1,PHI2)
MAC = ((PHI1'*PHI2).^2)./(diag(PHI1'*PHI1)*diag(PHI2'*PHI2)');
end
