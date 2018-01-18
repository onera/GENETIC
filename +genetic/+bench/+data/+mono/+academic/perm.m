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
function bench = perm(n)
if nargin == 0
   bench = [];
   return
end
% PERM
b        = 0.5;
bench    = struct('f'      , @(x) permFun(x, b, n),...
                  'bounds' , repmat([-n n],n,1));
end

function y = permFun(x, b, n)
s_out = 0;
for k=1:n
   s_in = 0;
   for j=1:n
      s_in = s_in + (j^k + b) * ( (x(j)/j)^k - 1);
   end
   s_out = s_out + s_in^2;
end
y = s_out;
end
