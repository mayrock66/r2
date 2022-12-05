function [p] = parser(fin)
%
% function [p] = parser(fin)
%
% September 2000
%
% Author: Jeffery L. Gray, grayj@purdue.edu
%         School of Electrical and Computer Engineering
%         Purdue University, West Lafayette, IN 47907
%__________________________
%
% Input -
%
%       fin: file id of input file containting statements
%
% Output -
%
%       p: a structure, defined below
%
% The following variables are always returned:
%
% p.ncard - is a string containing the name of the input card
% p.err - error messages, as follows:
%    normal return conditions:
%         p.err=-1; p.errmess='end of file';
%         p.err=1; p.errmess='no errors detected on assignment statement';
%         p.err=2; p.errmess='no errors detected on title statement';
%    error detected:
%         p.err=999; p.errmess=strcat('cannot mix numbers and strings 
%in: ',p.var(i).name);
%
% parser.m reads stements from the file fin which are of the form:
%
%       *titlecard  anything at all
%
%   or
%
%       cardname  var1=value1,var2=string2 var3=value3
%       +     array1=va1/va1/va3, array2=str1/str2/str3/str4
%       +     var4 value2
%
% Lines beginning with a blank or a $ are ignored.
%
% Variables have 3 types: number, string, or empty
%
% Statements may be any number of lines. The card name
% must start in column 1. The continuation symbol, +,
% must also appear in column 1.
%
% Commas or blanks are assumed to be separaters. Any
% number of separaters may appear between assignments.
% An assignment cannot contain any blanks, i.e.
%      block = - 12.0
% is not valid. It must be written as
%      block=-12.0   instead.
%
% Multiple values can be assigned to the variable
% by separating the values by /'s. For example
% wl=.3/.4/.5/.6/.7/.8
%    and
% ohmic=yes/no/no/no/yes
%
%
% Examples-
%
%       |<-- column 1
%
%####   newcard  x=3.0,rat=cat  alpha=1.0e12
%       +  wl=.3/.4/.5/inf/1e500 string=ab/cde/fghi/j purdue
%
%----   p.ncard: 'newcard'
%       p.nvar: 6
%       p.err: 1
%       p.errmess: 'no errors detected on assignment statement'
%       p.var: [1x6 struct]
%          p.var(1).name: 'x'
%          p.var(1).type: 'number'
%          p.var(1).nval: 1
%          p.var(1).val: 3
%          p.var(2).name: 'rat'
%          p.var(2).type: 'string'
%          p.var(2).nval: 1
%          p.var(2).val: {'cat'}
%          p.var(3).name: 'alpha'
%          p.var(3).type: 'number'
%          p.var(3).nval: 1
%          p.var(3).val: 1.0000e+012
%          p.var(4).name: 'wl'
%          p.var(4).type: 'number'
%          p.var(4).nval: 5
%          p.var(4)val: [0.3000 0.4000 0.5000 inf inf]
%          p.var(5).name: 'string'
%          p.var(5).type: 'string'
%          p.var(5).nval: 4
%          p.var(5).val: {'ab'  'cdcd'  'fgh'  'z'}
%          p.var(6).name: 'purdue'
%          p.var(6).type: 'empty'
%          p.var(6).nval: 0
%          p.var(6).val: []
%
%####   *mess  hello world!
%
%----   p.ncard: '*mess'
%       p.err: 2
%       p.errmess: 'no errors detected on title statement'
%       p.title: 'hello world!  '
%
%####   test2 f=false r=3/yes
%
%----   p.ncard: 'test2'
%       p.err: 999
%       p.errmess: 'cannot mix numbers and strings in: r'
%       p.nvar: 2
%       p.var: [1x2 struct]
%          p.var(1).name: 'f'
%          p.var(1).type: 'string'
%          p.var(1).nval: 1
%          p.var(1).val: {'false'}
%          p.var(2).name: 'r'
%          p.var(2).type: 'number'
%          p.var(2).nval: 2
%          p.var(2).val: 3

% Read first line of statement definition - ignore comment lines

line(1)=' ';

while (line(1) == ' ')

    if(feof(fin) == 1)
       p.ncard='End-of-File';
       p.err=-1;
       p.errmess='end of file';
       return;
    end

    line = fgetl(fin);
    if(line == -1)
       p.ncard='End-of-File';
       p.err=-1;
       p.errmess='end of file';
       return;
    end

    if(isempty(line)~=1)
       if(size(line)==[1 0] |line(1)=='$' | line(1)=='%')
          line(1)=' ';
       end
    else
       line(1)=' ';
    end

end

statement = line;

% Read remaining lines of definition statement if continued

line(1)=' ';
while (line(1) == ' ')

    if(feof(fin)==1)
       break;
    end

    line = fgetl(fin);

    if(size(line)==[1 0])
       line(1)=' ';
    end
  
    if(isempty(line)~=1)
      if(line(1)=='+')
         line(1)=' ';
         statement=strcat(statement,line);
      else
         fseek(fin,-(max(size(line))+1),'cof');
         break;
      end
    else
       break;
%      line(1)=' ';
    end
    

end

% replace all occurances of ',' with a space, ' '

statement = strrep(statement,',',' ');

% decode the statement

% check for title statement

if statement(1) == '*'
    [p.ncard rem]=strtok(statement);
    p.title=strjust(rem,'left');
    p.err=2;
    p.errmess='no errors detected on title statement';
    return;
end

% get statement name

[p.ncard rem]=strtok(statement);
rem=deblank(rem);

% If the end card is encounterd exit out from the parser

if strcmpi(p.ncard,'end')==1
   p.err=1;
   return;
end

% extract tokens

i=0;

while i>=0

    [token rem]=strtok(rem);
    n=size(token,2);
    if(n==0)
       break;
    end

    i=i+1;
    if i==1
       list=token;
    else
       list=char(list,token);
    end

end

% interpret list of tokens

n=size(list,1);
p.nvar=n;

for i=1:n

    [vname rem]=strtok(list(i,:),'=');
    rem(1,1)=' ';rem=strjust(rem,'left');rem=deblank(rem);
    p.var(i).name=deblank(vname);
    p.var(i).type='      ';
    if isempty(rem)==0
       % test for array of values
       delims=findstr(rem,'/');
       [dum nd]=size(delims);nv=nd+1;
       p.var(i).nval=nv;
       for j=1:nv
          [value rem]=strtok(rem,'/');
          rem(1,1)=' ';rem=strjust(rem,'left');rem=deblank(rem);
          number=str2double(value);
          if isnan(number)==1
             if p.var(i).type == 'number'
                p.err=999;
                p.errmess=strcat('cannot mix numbers and strings in:',p.var(i).name);
                return;
             end
             p.var(i).type='string';
             p.var(i).val{j}=value;
          else
             if p.var(i).type == 'string'
                p.err=999;
                p.errmess=strcat('cannot mix numbers and strings in:',p.var(i).name);
             end
             p.var(i).type='number';
             p.var(i).val(j)=number;
          end
       end
    else
       p.var(i).nval=0;
       p.var(i).type='empty';
    end

end

p.err=1;
p.errmess='no errors detected on assignment statement';

return;
