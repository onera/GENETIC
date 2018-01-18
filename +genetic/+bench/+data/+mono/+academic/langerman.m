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
function bench = langerman(n)
if nargin == 0
   bench = [];
   return
end
% LANGERMAN
a        = [9.681 0.667;9.4 2.041;8.025 9.152;2.196 0.415;8.074 8.777];
c        = [0.806 0.517 0.1 0.908 0.965];
bench    = struct('f'      , @(x) langermanFun(x, a, c, n),...
                  'bounds' , repmat([0 10],n,1));
end
function y = langermanFun(x, a, c, n)
y = 0;
for i=1:5
   s = 0;
   for j=1:n
      s = s + (x(j) - a(i,j))^2;
   end
   y = y - c(i)*exp(-s/pi)*cos(pi*s);
end
end
