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
function bench = lz09F6()
% lz09F6 benchmark
%
% Related reference
%  [1] H. Li and Q. Zhang, Multiobjective Optimization Problems with
%      Complicated Pareto Sets, MOEA/D and NSGA-II. IEEE Transactions
%      on Evolutionary Computation, vol. 2(12), pp. 284-302, 2009.
%

% Bench parameters
M = 3;
%
n = 10;
% Output
bench    = struct('f'      , @(x) lz09F6Fun(x),...
                  'fDim'   , M,...
                  'bounds' , [repmat([0 1],2,1);repmat([-2 2],n-2,1)],...
                  'fopt'   , lz09F6Front(M),...
                  'xDim'   , n);

end
% Handle for computation
function y = lz09F6Fun(x)
n           = length(x);
J1          = (4:3:n)';
J2          = (5:3:n)';
J3          = (3:3:n)';
y           = zeros(3,1);
y(1,1)      = cos(0.5*x(1,1)*pi)*cos(0.5*x(2,1)*pi) + 2*mean((x(J1,1) - 2*x(2,1)*ones(length(J1),1).*sin(2*pi*x(1,1)*ones(length(J1),1) + J1*pi/n)).^2,1);
y(2,1)      = cos(0.5*x(1,1)*pi)*sin(0.5*x(2,1)*pi) + 2*mean((x(J2,1) - 2*x(2,1)*ones(length(J2),1).*sin(2*pi*x(1,1)*ones(length(J2),1) + J2*pi/n)).^2,1);
y(3,1)      = sin(0.5*x(1,1)*pi) + 2*mean((x(J3,1) - 2*x(2,1)*ones(length(J3),1).*sin(2*pi*x(1,1)*ones(length(J3),1) + J3*pi/n)).^2,1);
end

function F = lz09F6Front(M)
F = genetic.bench.loadFrontData(['lz09F6Data' num2str(M) 'D'])';
end