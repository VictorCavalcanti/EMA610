function [D,JJ,DD,thetaf]=xupdate(Mm,Km,theta0,dMim,dKim)
%
%   Created by:     Daniel Kammer            
%                   Professor
%                   Dept of Engineering Physics
%                   University of Wisconsin
%                   Madison, WI  53706
%
%  Template of main function for finite element model updating.
%
%  Assumes mass and stiffness based design variables.
%
%
%
%  History:
%  =======
%
%
%=============================================================
%
%  n     = number of modal dof
%
%  ns    = number of sensors
%
%  na    = number of inputs
%
%  nvars = number of design variables
%
%
%  Input
%  =====
% 
%  Mm     = Modal mass matrix          (n x n)
%
%  Km     = Modal Stiffness matrix
%
%  f      = Frequency data matrix      (nd x 1)
%
%  theta0 = Nominal design variables   (nvars x 1)
%
%  dMim   = Mass sensitivity matrices normalized to dimensionless
%           design variables           (n x n x nvars)
%
%  dKim   = Stiffness sensitivity matrices normalized to dimensionless
%           design variables           (n x n x nvars)
%
%
%  Output
%  ======
%
%  thetaf  = final design variables   (nvars x 1)
%
%  DD      = fractional updated design variables by iteration
%
%  JJ      = cost function by iteration
%
%=================================================================
%
%
%  Usage:    [D,JJ,DD,thetaf]=xupdate(Mm,Km,theta0,dMim,dKim);
%

%  Maximum Number of Iterations
%  ----------------------------
fprintf('\n\n');
% maxnit = input('Maximum Number of Iterations - ') % prompt for maximum number of iterations
maxnit = 12;
alpha2=.01;                                   % design variable penalty number in J
pstop=.1;                                     % stopping criterion
n = length(diag(Mm));                         % number of modes
nvars = length(theta0);                       % number of design variables
nmodes = n;
ns = size(Mm,1);    
%  Set Weighting Matrices We, Wd
%  -----------------------------
We=eye(ns);                                   % sensor weighting matrix 
Wd=eye(nvars);                                % design variable weighting matrix

%  Set Bounds for Design Variable Changes During FMINCON Iterations
%  ----------------------------------------------------------------
DL=-.20*ones(nvars,1);
DU=.20*ones(nvars,1);

%**************************************
% ic = [1:na]';                                 % input column index
% in = [1:n]';                                  % mode number index

D = ones(nvars,1);                            % initial values of nondim design variables
nit = 1;                                      % counts number of iterations
DD = zeros(nvars,maxnit);                     % allocate space
JJ = zeros(maxnit,1);                         % allocate space

J = 10;
Jlast = 1;

%
%**************************************************************************
%**************************************************************************
    
while max(abs((J-Jlast)./Jlast))>pstop/100    % loop due to nonlinearity of sensitivities
    
    if nit>maxnit;                            % check for maximum iterations
        break
    end
    
    nit
    
    % Compute mode shapes/frequencies for current update
    
    [p,L] = eig(Km,Mm);                       % find mode shapes/frequencies of modal matrices
    MM = p'*Mm*p;
    p = p*diag(diag(MM).^-.5);                % mass-normalize modes
    [Ls,ix] = sort(diag(L));
    L = Ls;
    p = p(:,ix);                              % sort modes in order of increasing freq
    wa = sqrt(L);                             % natural frequencies, rad/s
    
    deps=zeros(ns,nvars,nf);                  % allocate space for depsilon/dD
   
            
    % Compute current sensitivities - depsilon/dD
    % pg 17 - https://ay16-17.moodle.wisc.edu/prod/pluginfile.php/171672/mod_resource/content/5/eCOWI_Resources/notes/lecture18-/Topic%2014%20-%20FEM%20Calibration%20Pres.pdf
    
    
    eps0 = (Pa-Pe)/b;                                                % compute error vector epsilon0_i
    
    Dlast = D;                                                       % value of D from previous iteration
    Jlast = J;                                                       % value of J from previous iteration
    
    [dD,J]=xmin(deps,We,Wd,eps0,alpha2,DL,DU);                       % find update dD
    
    D = D+dD                                                         % update D
    JJ(nit) = J;                                                     % cost function value for current iteration
    DD(:,nit) = D;                                                   % augment design variable history
    
    for k=1:nvars
        Km = Km+dKim(:,:,k)*dD(k);                                   % update modal stiffness matrix
        Mm = Mm+dMim(:,:,k)*dD(k);                                   % update modal mass matrix
    end
    
    thetaf = D.*theta0;                                              % final dimensional design vars
    
    nit = nit+1;                                                     % increment number of iterations

end

JJ = JJ(1:nit-1);
DD = DD(:,1:nit-1);



