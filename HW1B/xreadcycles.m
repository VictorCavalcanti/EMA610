function [cycles]=xreadcycles(dum)
%
%   Created by:     Ben Quasius
%                   Research Assistant
%                   Dept of Engineering Physics
%                   University of Wisconsin
%                   Madison, WI  53706
%
%  This function reads the cycles from the eigenvalue table print out found on a
%  NASTRAN .f06 file
%  ======================================================================
%
%  HISTORY
%  =======
%
%  Created:  Jun-26-00
%  Updated:  May-28-03
%
%  ======================================================================
%
%  INPUT
%  =====
%  dum       =  Dummy input
%
%  OUTPUT
%  ======
%  cycles	=  Vector containing cycles of modes
%               
%
%  REQUIRED FILES
%  ==============
%
%    NASTRAN .F06 file containing eigenvalue table
%
%
%  Use:  [cycles]=xreadcycles(dum);
%
%=======================================================================
%

f_in=0;
dof = zeros(1);
fprintf('\n')
filename2 = input('.f06 Filename Containing A-set Data:  ', 's'); % return as text variable

[f_in,message] = fopen(filename2,'r');

% opens the file 'filename2' 
if f_in == -1										 % checks success of fopen
   disp(message)
end

fseek(f_in,0,-1);									 % rewinds the file
%
% Find the location in file where 'U S E T' is first found
% --------------------------------------------------------
while 1
   mdata = fgetl(f_in);								 % data is a dummy variable
%   disp(mdata)
   match = findstr('NUMBER OF ROOTS FOUND',mdata);
   if length(mdata)<length('NUMBER OF ROOTS FOUND');             % makes sure mdata is longer than string
   elseif (isempty(match)==0);
      fseek(f_in, - length(mdata), 'cof');			     % rewinds to start of the line;
      [dum,a]=fscanf(f_in,'%s %s %s %s %s',5);            % scans current line for the 1st 5 strings
      [roots,b]=fscanf(f_in,'%i',1) ;                     % continues scan for one integer (# roots)
      break;
   end
end

cycles=0;
%
% Search the file for the column data
% -----------------------------------
k=1;       % k is a dummy variable that tracks the last mode read 
while 1
   if (feof(f_in) ~= 0)     % check end of file
      break;
   end
   % repeats till end of file or row is repeated
   [dum,a]=fscanf(f_in,'%i %i %g %g %g %g %g',[1 7]);
   if a<7                       % a<7 means that line wasn't a data line
      data = fgetl(f_in);       % reads remainder of line
   elseif a==7                  % a==7 successful line read.  
      mode=dum(1);              % first column is the mode number
      if mode<k                 % if mode number is less than k, a new set of data was reached
         break;                 % consequently, break out of loop
      end
      cycles(mode) = dum(5);    % assign the cycle value of the line to correct mode # in output vector
      k=mode;                   % set k to successfully read mode
      if mode==roots            % if mode==roots, then all modes have been read, so break
         break;
      end      
  end
end

cycles = cycles';
fclose(f_in);										 % Close the file
