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
function bench = zdt4()
n = 10;
lb    = -5 * ones(n,1);
ub    =  5 * ones(n,1);
lb(1) = 0;
ub(1) = 1;
% Output
bench    = struct('f'      , @(x) zdt4Fun(x),...
                  'fDim'   , 2,...
                  'bounds' , [lb,ub],...
                  'fopt'   , zdt4Front(),...
                  'xDim'   , n);

end
% handle for computation
function y = zdt4Fun(x)
y     = zeros(2,1);
y(1)  = x(1);
g     = 91 + sum(x(2:end).^2 - 10*cos(4*pi*x(2:end)));
h     = 1 - sqrt(y(1)/g);
y(2)  = g * h;
end

function F = zdt4Front()
% Obtained from jMetal
F = genetic.bench.loadFrontData('zdt4Data')';
end