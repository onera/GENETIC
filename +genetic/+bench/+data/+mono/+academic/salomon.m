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
function bench = salomon(n)
if nargin == 0
   bench = [];
   return
end
% SALOMON
a        = 2*pi;
b        = 0.1;
c        = 1;
bench    = struct('f'      , @(x)salomonFun(x,a,b,c),...
                  'bounds' , repmat([-4 4],n,1));
end

function y = salomonFun(x,a,b,c)
x = sqrt(x'*x);
y = -cos(a*nox(x)) + b*nox(x) + c;
end
