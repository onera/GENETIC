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
function bench = dropWave()
n = 2;
% Output
bench    = struct('f'      ,@(x) -(1+cos(12*norm(x)))/(2+0.5*norm(x)^2),...
                  'bounds' , repmat([-5.12 5.12],n,1),...
                  'xopt'   , zeros(n,1),...
                  'fopt'   , -1,...
                  'xDim'   , 2);

end
