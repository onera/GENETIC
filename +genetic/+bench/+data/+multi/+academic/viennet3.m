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
function bench = viennet3()
% Viennet3 benchmark
%

% Bench parameters
M = 3;
%
n = 2;
% Output
bench    = struct('f'      , @(x) viennet3Fun(x),...
                  'fDim'   , M,...
                  'bounds' , repmat([-3 3],n,1),...
                  'fopt'   , viennet3Front(M),...
                  'xDim'   , n);

end
% Handle for computation
function y = viennet3Fun(x)
y        = zeros(3,1);
y(1,1)   = 0.5*(x(1,1)^2 + x(2,1)^2) + sin(x(1,1)^2 + x(2,1)^2);
y(2,1)   = ((3*x(1,1) - 2*x(2,1) + 4)^2)/8 + ((x(1,1) - x(2,1) + 1)^2)/27 + 15;
y(3,1)   = 1/(x(1,1)^2 + x(2,1)^2 + 1) - 1.1*exp(-(x(1,1)^2 + x(2,1)^2));
end

function F = viennet3Front(M)
F = genetic.bench.loadFrontData(['viennet3Data' num2str(M) 'D'])';
end