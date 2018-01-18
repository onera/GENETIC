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
function bench = viennet2()
% Viennet2 benchmark
%

% Bench parameters
M = 3;
%
n = 2;
% Output
bench    = struct('f'      , @(x) viennet2Fun(x),...
                  'fDim'   , M,...
                  'bounds' , repmat([-4 4],n,1),...
                  'fopt'   , viennet2Front(M),...
                  'xDim'   , n);

end
% Handle for computation
function y = viennet2Fun(x)
y        = zeros(3,1);
y(1,1)   = ((x(1,1) - 2)^2)/2 + ((x(2,1) + 1)^2)/13 + 3;
y(2,1)   = ((x(1,1) + x(2,1) - 3)^2)/36 + ((-x(1,1) + x(2,1) + 2)^2)/8 - 17;
y(3,1)   = ((x(1,1) + 2*x(2,1) - 1)^2)/175 + ((2*x(2,1) - x(1,1))^2)/17 - 13;
end

function F = viennet2Front(M)
F = genetic.bench.loadFrontData(['viennet2Data' num2str(M) 'D'])';
end