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
function bench = cec09UF5()
% UF5 - CEC'2009 competition benchmark
%

% Bench parameters
M        = 2;
epsilon  = 0.1;
N        = 10;
%
n = 30;
% Output
bench    = struct('f'      , @(x) cec09UF5Fun(x, epsilon, N),...
                  'fDim'   , M,...
                  'bounds' , [[0 1];repmat([-1 1],n-1,1)],...
                  'fopt'   , cec09UF5Front(M),...
                  'xDim'   , n);

end
% Handle for computation
function y = cec09UF5Fun(x, epsilon, N)
n        = length(x);
J1       = (3:2:n)';
J2       = (2:2:n)';
tmp      = x - sin(6*pi*x(1,1)*ones(n,1) + (1:n)'*pi/n);
h        = 2*tmp.^2 - cos(4*pi*tmp) + 1;
y(1,1)   = x(1,1) + (1/(2*N) + epsilon)*abs(sin(2*N*pi*x(1,1))) + 2*mean(h(J1,1),1);
y(2,1)   = 1 - x(1,1) + (1/(2*N) + epsilon)*abs(sin(2*N*pi*x(1,1))) + 2*mean(h(J2,1),1);
end

function F = cec09UF5Front(M)
F = genetic.bench.loadFrontData(['cec09UF5Data' num2str(M) 'D'])';
end