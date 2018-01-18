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
function T = vonNeumann(n)
T     = cell(n,1);
igrid = floor(sqrt(n))+1;
jgrid = floor(n/igrid)+1;
for j = 1:n
   % Temporary vector <tmp> of dimension 5
   tmp = zeros(5,1);
   % From series to grid: <i> ---> (<ilig>,<jcol>)
   ilig = fix((j-1)/jgrid)+1;
   jcol = mod((j-1),jgrid)+1;
   % From grid to series: (<ilig>,<jcol>) ---> <i> = (<ilig>-1)*<jgrid> + <jcol>
   % >> [1ST] << %
   if jcol ~= jgrid
      tmp(1) = (ilig - 1)*jgrid + jcol + 1;
   else
      tmp(1) = (ilig - 1)*jgrid + 1;
   end
   % Check [1ST]
   if tmp(1) > n
      tmp(1) = (igrid - 1)*jgrid + 1;
   end
   % >> [2ND] << %
   if ilig ~= igrid
      tmp(2) = ilig*jgrid + jcol;
   else
      tmp(2) = jcol;
   end
   % Check [2ND]
   if tmp(2) > n
      tmp(2) = jcol;
   end
   % >> [3RD] << %
   if jcol ~= 1
      tmp(3) = (ilig - 1)*jgrid + jcol - 1;
   else
      tmp(3) = (ilig - 1)*jgrid + jgrid;
   end
   % Check [3RD]
   if tmp(3) > n
      tmp(3) = n;
   end
   % >> [4TH] << %
   if ilig ~= 1
      tmp(4) = (ilig - 2)*jgrid + jcol;
   else
      tmp(4) = (igrid - 1)*jgrid + jcol;
   end
   % Check [4TH]
   if tmp(4) > n
      tmp(4) = tmp(4) - jgrid;
   end
   % >> [5TH] << %
   tmp(5)   = j;
   % Update output idxList
   T{j}     = unique(setdiff(tmp,j));
end
end