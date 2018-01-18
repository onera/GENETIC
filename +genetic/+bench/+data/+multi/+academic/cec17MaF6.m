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
function bench = cec17MaF6()
% MaF6 - CEC'2017 competition benchmark
%
%  This test problem is used to assess whether evolutionary many-objectives
%  algorithms are capable of dealing with degenerate Pareto fronts.
%

% Bench parameters
M  = 3;
K  = 10;
I  = 2;
%
n = M-1 + K;
% Output
bench    = struct('f'      , @(x) cec17MaF6Fun(x, M, I),...
                  'fDim'   , M,...
                  'bounds' , repmat([0 1],n,1),...
                  'fopt'   , [],...
                  'xDim'   , n);

end
% Handle for computation
function y = cec17MaF6Fun(x, M, I)
g           = sum((x(M:end,1) - 0.5).^2,1);
tmp         = repmat(g,M-I,1);
x(I:M-1,1)  = (1 + 2*tmp.*x(I:M-1,1))./(2 + 2*tmp);
y           = flipud(cumprod([1;cos(x(1:M-1,1)*pi/2)],1)).*[1;sin(x(M-1:-1:1,1)*pi/2)].*repmat(1+100*g,M,1);
end

% function F = cec17MaF6Front(M)
% F = genetic.bench.multi.loadData(['cec17MaF6Data' num2str(M) 'D'])';
% end