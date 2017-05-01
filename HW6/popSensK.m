function dKim = popSensK
dKim = zeros(6,6,16);

%k1
dKim(:,:,1) = ...
 [1 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0];
%k2
dKim(:,:,2) = ...
 [1 -1 0 0 0 0;...
  -1 1 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0];
%k3
dKim(:,:,3) = ...
 [1 0 -1 0 0 0;...
  0 0 0 0 0 0;...
  -1 0 1 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0];
%k4
dKim(:,:,4) = ...
 [1 0 0 0 0 -1;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  -1 0 0 0 0 1];
%k5
dKim(:,:,5) = ...
 [0 0 0 0 0 0;...
  0 1 0 0 0 -1;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 -1 0 0 0 1];
%k6
dKim(:,:,6) = ...
 [0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 1 0 0 -1;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 -1 0 0 1];
%k7
dKim(:,:,7) = ...
 [0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 1 0 -1;...
  0 0 0 0 0 0;...
  0 0 0 -1 0 1];
%k8
dKim(:,:,8) = ...
 [0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 1 -1;...
  0 0 0 0 -1 1];
%k9
dKim(:,:,9) = ...
 [0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 1 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0];
%k10
dKim(:,:,10) = ...
 [0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 0 0;...
  0 0 0 0 1 0;...
  0 0 0 0 0 0];
end