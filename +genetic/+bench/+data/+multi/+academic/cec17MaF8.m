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
function bench = cec17MaF8()
% MaF8 - CEC'2017 competition benchmark
%
%  One important characteristic of 'cec17MaF8' is its Pareto optimal region
%  in the decision variables space which is typically a 2D manifold (regardless
%  of the dimensionality of its objective vectors). This naturally allows a
%  direct observation of the search behavior of evolutionary many-objectives
%  optimization algorithms, e.g., the convergence of their population to the
%  Pareto optimal solutions and the coverage of the population over the optimal
%  region.
%

% Bench parameters
M                    = 3;
[theta, rho]         = cart2pol(0,1);
[pts(:,1), pts(:,2)] = pol2cart(theta - 2*pi*(1:M)'/M, rho);
pts                  = pts';
%
n = 2;
% if n ~= 2
%    error('The benchmark ''cec17MaF8'' is only available for dimension 2.')
% end
% Output
bench    = struct('f'      , @(x) cec17MaF8Fun(x, pts),...
                  'fDim'   , M,...
                  'bounds' , repmat([-1e4 1e4],n,1),...
                  'fopt'   , [],...
                  'xDim'   , 2);

end
% Handle for computation
function y = cec17MaF8Fun(x, pts)
y = genetic.tools.pairwiseDistance(pts, x);
end

% function F = cec17MaF8Front(M)
% F = genetic.bench.multi.loadData(['cec17MaF8Data' num2str(M) 'D'])';
% end