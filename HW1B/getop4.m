function [name, M]=getop4(f_id)
%
% This function takes the already open file specified by f_id 
% (f_id is an integer file identifier obtained from "fopen")
% and finds the next op4 matrix in the file.  It returns the matrix 
% and the name of the matrix.  This function can then be called 
% repeatedly to read successive matrices
% ===================================================================
%
%  HISTORY
%  =======
%
%  Created:  Mar-20-00    Ben Quasius
%  Modified: Jun-07-00    Ben Quasius     Bug fix
%
%  =======================================================================
%
%  INPUT
%  =====
%  f_id      =  Interger file identifier of .OP4 file obtained from "fopen"
%
%  OUTPUT
%  ======
%  name      =  Name of matrix read
%
%  M         =  Matrix read
%
%
%  Use:  [name, M]=getop4(f_id);
%
%==================================================================
%
%fseek(f_id,0,-1);					        % rewinds the file
%
position=ftell(f_id);
head = fscanf(f_id,'%i %i %i %i',[4 1]);   	% [N, M, I2, I2] %
nam = fscanf(f_id, '%s32')				    % matrix name
position2=ftell(f_id);
fseek(f_id, position-position2,0);
fgetl(f_id);
%
% Next we must import the matrix.  The size head(2) by head(1)
% ------------------------------------------------------------
A=zeros(head(2),head(1));	 % zeroes matrix MxN
%
% repeat scanning file until col num is larger than N 
% ---------------------------------------------------
%col=0;							            % initalize col index
%
while(1==1);					
   colhead=fscanf(f_id,'%i %i %i',[3 1]);	% read col headers
   col=colhead(1);
   if (col>head(1)); 
      break; 
   end;      
   for i=0:colhead(3)-1;				      % populate A matrix
      B=fscanf(f_id,'%e',[1,1]);
      A(colhead(2)+i,colhead(1))=B(1,1);
   end
end
name = nam;
M = A;
size(M)




 
