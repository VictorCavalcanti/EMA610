function [dof]=xreadasetA(dum)
%
%   Created by:     Ben Quasius
%                   Research Assistant
%                   Dept of Engineering Physics
%                   University of Wisconsin
%                   Madison, WI  53706
%
%  This function reads the A-set from the USET table print out found on a
%  NASTRAN .f06 file for NASTRAN Version 70.5
%  ======================================================================
%
%  HISTORY
%  =======
%
%  Created:  Apr-05-00
%  Modified: Jun-14-00     Bug fix
%  Modified  Oct-19-00     Read A set instead of BF set.    DCK
%  Modified  Oct-15-01     Bug fixes AND documentation.    BEQ
%  Modified  Mar-15-02     Bug fixes BEQ
%  ======================================================================
%
%  INPUT
%  =====
%  dum       =  Dummy input
%
%  OUTPUT
%  ======
%  dof       =  Vector containing DOF ID numbers in numerical order
%               from A-set
%               ie.  8000-1 = 8000.1 etc.
%
%
%  REQUIRED FILES
%  ==============
%
%    NASTRAN .F06 file containing USET table
%
%
%  Use:  [dof]=xreadasetA(dum);
%
%  NOTES
%  =====
%
%  data			dummy string of text from file
%  k, row & q	dummy index variables
%  match		checks 'data' strings for line data
%
%  BUG FIXES
%  =========
%  10/12/01  BEQ
%  Fixed bug that adds an extra degree of freedom for some titles.  
%  Case in point: The Nastran job title "X33 ..." would 
%  add 4.3 as the final dof. This was fixed by checking that the 
%  character skipped was a "-".  Any other character will break
%  the fscan while loop.  The program will continue to search for 
%  lines that contain "#1=" and add dofs until the number is 
%  smaller than the dof set size (signifying a new data set).
%
%  Also, the data I/O has been streamlined to read and rewind 
%  less frequently.  This is accomplished by using sscabf instead
%  of fscanf and fseek.
%
%  3/13/02  BEQ
%  Fixed bug that reads the first set of data.  The script now only 
%  checks one line for 'A       DISPLACEMENT'.  
%  Also, "space ~= '-'" has been changed to "not(strcmp(space,'-'))" 
%  to suppress warning.
%=======================================================================
%
% Open the F06 file
f_in=0;
dof = zeros(1);
fprintf('\n')
filename2 = input('.f06 Filename Containing A-set Data:  ', 's'); % return as text variable
[f_in,message] = fopen(filename2, 'r');		         % opens the file 'filename2'
if f_in == -1		 % checks success of fopen
    disp(message)
end
fseek(f_in,0,-1);		 % rewinds the file
%
% Find the location in file where 'U S E T' is first found
% -----mda  ---------------------------------------------------
mdata = ' '; % initializing mdata 
%
%  Search for 'A    DISPLACEMENT SET'
%  If the f06 file is word wrapped (on propose or accident)
%  'A' and 'DISPLACEMENT SET' could be on different lines
%  This script directly compares the strings, checks two
%  lines at a time, and removes spaces in the strings and checks
%  again.
while 1;  
    mdata = fgetl(f_in);   %get new line
    test=sscanf(mdata,'%s');   %strip spaces from line
    match = strcmp('ADISPLACEMENTSET',test);
    if (length(mdata)<length('A        DISPLACEMENT SET'));                   
        %makes sure mdata is longer than string
    elseif match==1;
        break;
    end
end
%
% Search the file for the column data
% -----------------------------------
k=1; row=1;  %initialize k and row
%  This loop terminates when the end of file is reached or the row # is less 
%  than the number of rows in the DOF, i.e., a new set of data is reached.
while (row >= k)&(feof(f_in) == 0)	;				 
    % repeats till end of file or row is repeated
    data = fgetl(f_in);
    match = findstr(data, '1='); % check line for row indicator
    if match > 0;  %if the match is positive, start retreiving the dof
        fseek(f_in, -length(data)-2, 'cof'); %rewinds to start of the line
        row = fscanf(f_in,'%d',1); %extracts row number
        if row<k;
            break;
        end;
        eq_str=fscanf(f_in,'%s',1); % skips '=' string
        while row >= k  % analyzing each line
            [nod,a]=fscanf(f_in,'%d',1);
            [space,b]=fscanf(f_in,'%1s',1);
            if not(strcmp(space,'-'));break;end   %verify that skipped character is a dash
            [dir,c]=fscanf(f_in,'%d',1);
            %continue if node,dash,dir have been correctly read
            if ((a==0)+(b==0)+(c==0))>0  
                break; 
            end
            dof(row)=nod+dir/10;
            row=row+1;
        end	% end of inner while loop
        k=row;
        
    end	    % end of if
end         % end of while
%
dof = dof';
%fclose(f_in); % Close the file
