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
function bench = cec17MaF7()
% MaF7 - CEC'2017 competition benchmark
%
%  This test problem is used to assess whether evolutionary many-objectives
%  algorithms are capable of dealing with disconnected Pareto fronts, when
%  moreover the number of disconnected segments is large in high-dimensional
%  objective spaces.
%

% Bench parameters
M  = 3;
K  = 10;
%
n = M-1 + K;
% Output
bench    = struct('f'      , @(x) cec17MaF7Fun(x, M),...
                  'fDim'   , M,...
                  'bounds' , repmat([0 1],n,1),...
                  'fopt'   , [],...
                  'xDim'   , n);

end
% Handle for computation
function y = cec17MaF7Fun(x, M)
y           = zeros(M,1);
g           = 1 + 9*mean(x(M:end,1),1);
y(1:M-1,1)  = x(1:M-1,1);
y(M,1)      = (M - sum(x(1:M-1,1)./(1 + repmat(g,M-1,1)).*(1 + sin(3*pi*x(1:M-1,1))),1))*(1 + g);
end

% function F = cec17MaF7Front(M)
% F = genetic.bench.multi.loadData(['cec17MaF7Data' num2str(M) 'D'])';
% end