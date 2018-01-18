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
function bench = cec17MaF4()
% MaF4 - CEC'2017 competition benchmark
%
%  This test problem is used to assess whether evolutionary many-objectives
%  algorithms are capable of dealing with badly-scaled Pareto fronts notably
%  when the fitness landscape is highly multidmodal.
%

% Bench parameters
M = 3;
K = 10;
%
n = M-1 + K;
% Output
bench    = struct('f'      , @(x) cec17MaF4Fun(x, M, K),...
                  'fDim'   , M,...
                  'bounds' , repmat([0 1],n,1),...
                  'fopt'   , [],...
                  'xDim'   , n);

end
% Handle for computation
function y = cec17MaF4Fun(x, M, K)
g  = 100*(K + sum((x(M:end,1) - 0.5).^2 - cos(20*pi*(x(M:end,1) - 0.5)),1));
y  = repmat(1+g,M,1) - flipud(cumprod([1;cos(x(1:M-1,1)*pi/2)],1)).*[1;sin(x(M-1:-1:1,1)*pi/2)].*repmat(1+g,M,1);
y  = (2.^(1:M)').*y;
end

% function F = cec17MaF4Front(M)
% F = genetic.bench.multi.loadData(['cec17MaF4Data' num2str(M) 'D'])';
% end