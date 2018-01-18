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
function bench = cec09UF4()
% UF4 - CEC'2009 competition benchmark
%

% Bench parameters
M = 2;
%
n = 30;
% Output
bench    = struct('f'      , @(x) cec09UF4Fun(x),...
                  'fDim'   , M,...
                  'bounds' , [[0 1];repmat([-2 2],n-1,1)],...
                  'fopt'   , cec09UF4Front(M),...
                  'xDim'   , n);

end
% Handle for computation
function y = cec09UF4Fun(x)
n        = length(x);
J1       = (3:2:n)';
J2       = (2:2:n)';
tmp      = x - sin(6*pi*x(1,1)*ones(n,1) + (1:n)'*pi/n);
h        = abs(tmp)./(1 + exp(2*abs(tmp)));
y(1,1)   = x(1,1) + 2*mean(h(J1,1),1);
y(2,1)   = 1 - x(1,1)^2 + 2*mean(h(J2,1),1);
end

function F = cec09UF4Front(M)
F = genetic.bench.loadFrontData(['cec09UF4Data' num2str(M) 'D'])';
end