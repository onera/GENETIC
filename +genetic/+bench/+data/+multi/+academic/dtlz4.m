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
function bench = dtlz4()
M     = 2;
K     = 10;
alpha = 100;
n = M-1 + K;
% if n~=((M-1)+K)
%    error('The benchmark ''dtlz4'' is only available for dimension %d.',(M-1)+K)
% end
% Output
bench    = struct('f'      , @(x) dtlz4Fun(x, M, K, alpha),...
                  'fDim'   , M,...
                  'bounds' , repmat([0 1],n,1),...
                  'fopt'   , dtlz4Front(M),...
                  'xDim'   , n);

end
% handle for computation
function y = dtlz4Fun(x, M, K, alpha)
n           = (M - 1) + K;
xm          = x(n-K+1:end,1);
g           = sum((xm - 0.5).^2, 1);
y(1,1)      = (1 + g).*prod(cos(pi/2*x(1:M-1,:).^alpha),1);
for i = 2:M-1
   y(i,1) = (1 + g).*prod(cos(pi/2*x(1:M-i,:).^alpha),1).*sin(pi/2*x(M-i+1,:).^alpha);
end
y(M,1)      = (1 + g).*sin(pi/2*x(1,1).^alpha);
end

function F = dtlz4Front(M)
% Obtained from jMetal
F = genetic.bench.loadFrontData(['dtlz4Data' num2str(M) 'D'])';
end