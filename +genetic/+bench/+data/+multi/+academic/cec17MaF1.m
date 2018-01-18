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
function bench = cec17MaF1()
% MaF1 - CEC'2017 competition benchmark
%
%  This test problem is used to assess whether evolutionary many-objectives
%  algorithms are capable of dealing with inverted Pareto fronts.
%

% Bench parameters
M = 3;
K = 10;
%
n = M-1 + K;

% Output
bench    = struct('f'      , @(x) cec17MaF1Fun(x, M),...
                  'fDim'   , M,...
                  'bounds' , repmat([0 1],n,1),...
                  'fopt'   , [],...
                  'xDim'   , n);

end
% Handle for computation
function y = cec17MaF1Fun(x, M)
g        = sum((x(M:end,1) - 0.5).^2,1);
y(1,1)   = (1 - prod(x(1:M-1,1),1))*(1 + g);
for i = 2:M-1
   y(i,1) = (1 + g)*(1 - prod(x(1:M-i,1),1)*(1 - x(M-i+1,1)));
end
y(M,1)   = x(1,1)*(1 + g);
end

% function F = cec17MaF1Front(M)
% F = genetic.bench.multi.loadData(['cec17MaF1Data' num2str(M) 'D'])';
% end