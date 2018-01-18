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
function bench = lz09F7()
% lz09F7 benchmark
%
% Related reference
%  [1] H. Li and Q. Zhang, Multiobjective Optimization Problems with
%      Complicated Pareto Sets, MOEA/D and NSGA-II. IEEE Transactions
%      on Evolutionary Computation, vol. 2(12), pp. 284-302, 2009.
%

% Bench parameters
M = 2;
%
n = 10;
% Output
bench    = struct('f'      , @(x) lz09F7Fun(x),...
                  'fDim'   , M,...
                  'bounds' , repmat([0 1],n,1),...
                  'fopt'   , lz09F7Front(M),...
                  'xDim'   , n);

end
% Handle for computation
function y = lz09F7Fun(x)
n           = length(x);
J1          = (3:2:n)';
J2          = (2:2:n)';
y           = zeros(2,1);
expVector   = 1 + 3*((1:n)'-2)/(n-2)/2;
tmp         = x - x(1,1)*ones(n,1).^expVector;
y(1,1)      = x(1,1) + 2*mean(4*tmp(J1,1).^2 - cos(8*tmp(J1,1)*pi) + 1,1);
y(2,1)      = 1 - sqrt(x(1,1)) + 2*mean(4*tmp(J2,1).^2 - cos(8*tmp(J2,1)*pi) + 1,1);
end

function F = lz09F7Front(M)
F = genetic.bench.loadFrontData(['lz09F7Data' num2str(M) 'D'])';
end