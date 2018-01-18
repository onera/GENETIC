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
function bench = shaffer1()
n = 1;
% Output
bench    = struct('f'      , @(x) shaffer1Fun(x),...
                  'fDim'   , 2,...
                  'bounds' , repmat([-1e5 1e5],n,1),...
                  'fopt'   , shaffer1Front(),...
                  'xDim'   , n);

end
% handle for computation
function y = shaffer1Fun(x)
y     = zeros(2,1);
y(1)  = x^2;
y(2)  = (x-2)^2;
end

function F = shaffer1Front()
% Obtained from jMetal
F = genetic.bench.loadFrontData('shaffer1Data')';
end