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
function bench = cec09UF9(n)
% UF9 - CEC'2009 competition benchmark
%

% Bench parameters
M = 3;
%
n = 30;
% Output
bench    = struct('f'      , @(x) cec09UF9Fun(x),...
                  'fDim'   , M,...
                  'bounds' , [repmat([0 1],2,1);repmat([-2 2],n-2,1)],...
                  'fopt'   , cec09UF9Front(M),...
                  'xDim'   , n);

end
% Handle for computation
function y = cec09UF9Fun(x)
n        = length(x);
J1       = (4:3:n)';
J2       = (5:3:n)';
J3       = (3:3:n)';
tmp      = x - 2*x(2,1)*ones(n,1).*sin(2*pi*x(1,1)*ones(n,1) + (1:n)'*pi/n);
y(1,1)   = 0.5*(max(0,1.1*(1 - 4*(2*x(1,1) - 1)^2)) + 2*x(1,1))*x(2,1) + 2*mean(tmp(J1,1).^2,1);
y(2,1)   = 0.5*(max(0,1.1*(1 - 4*(2*x(1,1) - 1)^2)) - 2*x(1,1) + 2)*x(2,1) + 2*mean(tmp(J2,1).^2,1);
y(3,1)   = 1 - x(2,1) + 2*mean(tmp(J3,1).^2,1);
end

function F = cec09UF9Front(M)
F = genetic.bench.loadFrontData(['cec09UF9Data' num2str(M) 'D'])';
end