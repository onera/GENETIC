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
function bench = whitley(n)
if nargin == 0
   bench = [];
   return
end
% WHITLEY
a        = 100;
b        = 4000;
c        = 1;
bench    = struct('f'      , @(x) whitleyFun(x, a, b, c, n),...
                  'bounds' , repmat([-10 10],n,1));
end
function y = whitleyFun(x, a, b, c, n)
y = 0;
for i=1:n
   for j=1:n
      y = y + (a*(x(i)^2 - x(j))^2 + (1 - x(j))^2)^2/b - cos(a*(x(i)^2 - x(j))^2 + (1 - x(j))^2) + c;
   end
end
end
