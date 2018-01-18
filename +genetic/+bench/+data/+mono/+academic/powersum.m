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
function bench = powersum(n)
if nargin == 0
   bench = [];
   return
end
% POWER SUM
b        = [8 18 44 114];
bench    = struct('f'      , @(x) powersumFun(x, b, n),...
                  'bounds' , repmat([0 n],n,1));

end
% POWER SUM
function y = powersumFun(x, b, n)
s_out = 0;
for k=1:n
   s_in = 0;
   for j=1:n
      s_in = s_in + x(j)^k;
   end
   s_out = s_out + (s_in - b(k))^2;
end
y = s_out;
end
