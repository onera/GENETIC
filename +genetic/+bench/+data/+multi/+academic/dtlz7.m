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
function bench = dtlz7()
M = 2;
K = 20;
n = M -1 + K;
% if n~=((M-1)+K)
%    error('The benchmark ''dtlz7'' is only available for dimension %d.',(M-1)+K)
% end
% Output
bench    = struct('f'      , @(x) dtlz7Fun(x, M, K),...
                  'fDim'   , M,...
                  'bounds' , repmat([0 1],n,1),...
                  'fopt'   , dtlz7Front(M),...
                  'xDim'   , n);

end
% handle for computation
function y = dtlz7Fun(x, M, K)
n           = (M - 1) + K;
xm          = x(n-K+1:end,1);
g           = 1 + 9/K*sum(xm,1);
y(1:M-1,1)  = x(1:M-1,1);
gaux        = g(ones(M-1,1),1);
h           = M - sum(y./(1 + gaux).*(1 + sin(3*pi*y)),1);
y(M,1)      = (1 + g).*h;
end

function F = dtlz7Front(M)
% Obtained from jMetal
F = genetic.bench.loadFrontData(['dtlz7Data' num2str(M) 'D'])';
end