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
function bench = cec17MaF5()
% MaF5 - CEC'2017 competition benchmark
%
%  This test problem is used to assess whether evolutionary many-objectives
%  algorithms are capable of dealing with badly-scaled Pareto fronts/shapes.
%

% Bench parameters
M     = 3;
K     = 10;
alpha = 100;
%
n = M-1 + K;
% Output
bench    = struct('f'      , @(x) cec17MaF5Fun(x, M, alpha),...
                  'fDim'   , M,...
                  'bounds' , repmat([0 1],n,1),...
                  'fopt'   , [],...
                  'xDim'   , n);

end
% Handle for computation
function y = cec17MaF5Fun(x, M, alpha)
x(1:M-1,1)  = x(1:M-1,1).^alpha;
g           = sum((x(M:end,1) - 0.5).^2,1);
y           = flipud(cumprod([1;cos(x(1:M-1,1)*pi/2)],1)).*[1;sin(x(M-1:-1:1,1)*pi/2)].*repmat(1+g,M,1);
y           = (2.^(M:-1:1)').*y;
end

% function F = cec17MaF5Front(M)
% F = genetic.bench.multi.loadData(['cec17MaF5Data' num2str(M) 'D'])';
% end