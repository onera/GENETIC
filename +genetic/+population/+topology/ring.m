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
function T = ring(n)
% Returns a ring topology, i.e. a topology in which each individual is
% connected to the previous one and the next one.

T = cell(n,1);
for i = 1:n
   if i > 1
      previous = i-1;
   else
      previous = n;
   end
   if i < n
      next = i+1;
   else
      next = 1;
   end
   T{i}     = [previous;next];
end

end