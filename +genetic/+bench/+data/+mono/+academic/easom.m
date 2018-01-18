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
function bench = easom()
n = 2;
bench = struct('f'      , @(x) -cos(x(1))*cos(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2),...
               'bounds' , repmat([-100 100],n,1),...
               'fopt'   , -1,...
               'xopt'   , pi*ones(n,1),...
               'xDim'   , 2);
end
