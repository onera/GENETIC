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
function bench = fonseca(n)
% Output
if nargin == 0
   bench = [];
   return
end
bench    = struct('f'      , @(x) fonsecaFun(x,n),...
                  'fDim'   , 2,...
                  'bounds' , repmat([-4 4],n,1),...
                  'fopt'   , fonsecaFront());

end
% handle for computation
function y = fonsecaFun(x,n)
y     = zeros(2,1);
coeff = 1/sqrt(n);
y(1)  = 1 - exp(-sum((x - coeff).^2));
y(2)  = 1 - exp(-sum((x + coeff).^2));
end

function F = fonsecaFront()
% Obtained from jMetal
F = genetic.bench.loadFrontData('fonsecaData')';
end