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
function bench = cec09UF2()
% UF2 - CEC'2009 competition benchmark
%

% Bench parameters
M = 2;
%
n = 30;
% Output
bench    = struct('f'      , @(x) cec09UF2Fun(x),...
                  'fDim'   , M,...
                  'bounds' , [[0 1];repmat([-1 1],n-1,1)],...
                  'fopt'   , cec09UF2Front(M),...
                  'xDim'   , n);

end
% Handle for computation
function y = cec09UF2Fun(x)
n           = length(x);
J1          = (3:2:n)';
J2          = (2:2:n)';
tmp         = zeros(n,1);
x1          = x(1,1)*ones(length(J1),1);
tmp(J1,1)   = x(J1,1) - (0.3*(x1.^2).*cos(24*pi*x1 + 4*J1*pi/n) + 0.6*x1).*cos(6*pi*x1 + J1*pi/n);
x1          = x(1,1)*ones(length(J2),1);
tmp(J2,1)   = x(J2,1) - (0.3*(x1.^2).*cos(24*pi*x1 + 4*J2*pi/n) + 0.6*x1).*sin(6*pi*x1 + J2*pi/n);
y(1,1)      = x(1,1) + 2*mean(tmp(J1,1).^2,1);
y(2,1)      = 1 - sqrt(x(1,1)) + 2*mean(tmp(J2,1).^2,1);
end

function F = cec09UF2Front(M)
F = genetic.bench.loadFrontData(['cec09UF2Data' num2str(M) 'D'])';
end