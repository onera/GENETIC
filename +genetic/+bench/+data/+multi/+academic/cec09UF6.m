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
function bench = cec09UF6()
% UF6 - CEC'2009 competition benchmark
%

% Bench parameters
M        = 2;
epsilon  = 0.1;
N        = 2;
%
n = 30;
% Output
bench    = struct('f'      , @(x) cec09UF6Fun(x, epsilon, N),...
                  'fDim'   , M,...
                  'bounds' , [[0 1];repmat([-1 1],n-1,1)],...
                  'fopt'   , cec09UF6Front(M),...
                  'xDim'   , n);

end
% Handle for computation
function y = cec09UF6Fun(x, epsilon, N)
n        = length(x);
J1       = (3:2:n)';
J2       = (2:2:n)';
tmp      = x - sin(6*pi*x(1,1)*ones(n,1) + (1:n)'*pi/n);
y(1,1)   = x(1,1) + max(0,2*(1/(2*N) + epsilon)*sin(2*N*pi*x(1,1))) + 2*(4*sum(tmp(J1,1).^2,1) - 2*prod(cos(20*tmp(J1,1)*pi./sqrt(J1)),1) + 2)/length(J1);
y(2,1)   = 1 - x(1,1) + max(0,2*(1/(2*N) + epsilon)*sin(2*N*pi*x(1,1))) + 2*(4*sum(tmp(J2,1).^2,1) - 2*prod(cos(20*tmp(J2,1)*pi./sqrt(J2)),1) + 2)/length(J2);
end

function F = cec09UF6Front(M)
F = genetic.bench.loadFrontData(['cec09UF6Data' num2str(M) 'D'])';
end