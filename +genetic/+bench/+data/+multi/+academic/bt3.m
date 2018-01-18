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
function bench = bt3()
% BT3 benchmark - Multi-objective test problems with bias
%
%  This multi-objective test problem  is characterized by two criteria that
%  must be optimized. The number of decision variables is set to 30.
%

% Bench parameters
n = 30;
M = 2;

% Output
bench    = struct('f'      , @(x) bt3Fun(x, M),...
                  'fDim'   , M,...
                  'bounds' , repmat([0 1],n,1),...
                  'fopt'   , [],...
                  'xDim'   , n);

end
% Handle for computation
function y = bt3Fun(x, M)
y     = zeros(M,1);
n     = length(x);
ind1  = (2:2:n);
ind2  = (3:2:n);
tmp   = x - sin(pi*(1:n)'/2/n);
y(1)  = abs(x(1,1))^0.02 + sum(tmp(ind1,1).^2 + (1 - exp(-tmp(ind1,1).^2/1e-8))/5,1);
y(2)  = 1 - sqrt(abs(x(1,1))^0.02) + sum(tmp(ind2,1).^2 + (1 - exp(-tmp(ind2,1).^2/1e-8))/5,1);
end

% function F = bt3Front(M)
% F = genetic.bench.multi.loadData(['bt3Data' num2str(M) 'D'])';
% end