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
function bench = rosenbrock(n)
if nargin == 0
   bench = [];
   return
end
% ROSENBROCK
a        = 100;
bench    = struct('f'      , @(x) sum(a*(x(1:n-1).^2 - x(2:n)).^2 + (x(1:n-1) - 1).^2),...
                  'bounds' , repmat([-5 10],n,1),...
                  'fopt'   , 0,...
                  'xopt'   , ones(n,1));
end
