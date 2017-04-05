function [name,M]=readop4(dum)
%
%   Created by:     Daniel C. Kammer
%                   Assistant Professor
%                   Dept of Engineering Mechanics
%                   University of Wisconsin
%                   Madison, WI  53706
%                   (608) 262-5724 / 262-3990
%
%  This function reads in a matrix from a nastran output4 BCD file
%  =================================================================
%
%  HISTORY
%  =======
%
%  Version 1.0  -  Created:  Dec-03-00   Daniel C. Kammer
%
%  INPUT
%  =====
%  dum       =  Dummy input set to anything
%
%  OUTPUT
%  ======
%  name      =  Matrix name on file
%
%  M         =  Matrix
%
%  --------------------------------------------------------
%
%  Use:    [name,M]=readop4(dum);
%
%  ========================================================
%
fprintf('\n\n')
filenameP = input('Matrix OP4 Filename:  ', 's');             % return as text variable
		[fid2,message] = fopen(filenameP, 'r');
		[name,M]=getop4(fid2);                                % read in mode matrix
		fprintf('\n')
%
