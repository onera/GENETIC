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
function bench = crossLegTable()
n        = 2;
fun      = @(x) -1/((abs(exp(abs(100-sqrt(x(1)^2 + x(2)^2)/pi)) * sin(x(1)) * sin(x(2)) ) + 1 )^0.1);
bench    = struct('f'      ,fun,...
                  'bounds' ,repmat([-10 10],n,1),...
                  'xopt'   , zeros(n,1),...
                  'fopt'   , -1,...
                  'xDim'   , n);
end

