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
function bench = zdt6()
n = 10;
% Output
bench    = struct('f'      , @(x) zdt6Fun(x),...
                  'fDim'   , 2,...
                  'bounds' , repmat([0 1],n,1),...
                  'fopt'   , zdt6Front(),...
                  'xDim'   , n);

end
% handle for computation
function y = zdt6Fun(x)
y     = zeros(2,1);
y(1)  = 1-exp(-4*x(1)) * sin(6*pi*x(1))^6;
g     = 1+9*sum((x(2:end)/9).^(0.25));
h     = 1 - (y(1)/g)^2;
y(2)  = g * h;
end

function F = zdt6Front()
% Obtained from jMetal
F = genetic.bench.loadFrontData('zdt6Data')';
end