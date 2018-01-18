% Copyright 2018 ONERA
%
% This file is part of the GENETIC project.
%
% GENETIC is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License version 3 as
% published by the Free Software Foundation.
%
% GENETIC is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with GENETIC.  If not, see <https://www.gnu.org/licenses/lgpl-3.0>.
%
function strOut = printCol(str, a, len)
ncol = length(str);
if nargin < 2 || isempty(a)
   a = repmat('r',1,ncol);
end
if nargin < 3
   len = 10*ones(ncol,1);
end
strOut = '';
margin = ' ';
for i = 1:ncol
   if i<ncol
      if isempty(str{i})
        str{i} = '';
      end
      strOut = [strOut, strAlign(a(i), str{i}, len(i)), margin];
   else
      strOut = [strOut, strAlign(a(i), str{i}, len(i))];
   end
end


end

function strOut = strAlign(side, str, n)
strLength = length(str);
if strLength > n
   str = str(1:n);
end
switch side
   case 'l'
      blankr      = printBlanks(n-strLength);
      strOut      = sprintf('%s%s',str, blankr);
   case 'r'
      blankl      = printBlanks(n-strLength);
      strOut      = sprintf('%s%s',blankl, str);
   case 'c'
      if mod(strLength,2)==0
         blankl  = printBlanks(round(n/2-strLength/2));
         blankr  = printBlanks(round(n/2-strLength/2));
         strOut  = sprintf('%s%s%s',blankl, str, blankr);
      else
         blankl  = printBlanks(round(n/2-(strLength-1)/2));
         blankr  = printBlanks(round(n/2-(strLength-1)/2-1));
         strOut  = sprintf('%s%s%s',blankl, str, blankr);
      end   
end

end

function strOut = printBlanks(n)
strOut = '';
if n > 0
   strOut   = blanks(n);
end
end