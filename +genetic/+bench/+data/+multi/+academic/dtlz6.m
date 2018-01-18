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
function bench = dtlz6()
M     = 2;
K     = 10;
n     = M -1 + K;
% if n~=((M-1)+K)
%    error('The benchmark ''dtlz6'' is only available for dimension %d.',(M-1)+K)
% end
% Output
bench    = struct('f'      , @(x) dtlz6Fun(x, M, K),...
                  'fDim'   , M,...
                  'bounds' , repmat([0 1],n,1),...
                  'fopt'   , dtlz6Front(M),...
                  'xDim'   , n);

end
% handle for computation
function y = dtlz6Fun(x, M, K)
n              = (M - 1) + K;
xm             = x(n-K+1:end,1);
g              = sum(xm.^0.1,1); 
theta(1,1)     = pi/2*x(1,1);
gr             = g(ones(M-2,1),:);
theta(2:M-1,1) = pi./(4*(1+gr)).*(1 + 2*gr.*x(2:M-1,1));
y(1,1)         = (1 + g).*prod(cos(theta(1:M-1,1)),1);
for i = 2:M-1
   y(i,1) = (1 + g).*prod(cos(theta(1:M-i,1)),1).*sin(theta(M-i+1,1));
end
y(M,1)         = (1 + g).*sin(theta(1,1));
end

function F = dtlz6Front(M)
% Obtained from jMetal
F = genetic.bench.loadFrontData(['dtlz6Data' num2str(M) 'D'])';
end