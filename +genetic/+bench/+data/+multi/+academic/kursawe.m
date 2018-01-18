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
function bench = kursawe()
n = 3;
% Output
bench    = struct('f'      , @(x) kursaweFun(x),...
                  'fDim'   , 2,...
                  'bounds' , repmat([-5 5],n,1),...
                  'fopt'   , kursaweFront(),...
                  'xDim'   , n);

end
% handle for computation
function y = kursaweFun(x)
y     = zeros(2,1);
y(1)  = sum(-10*exp(-0.2*(sqrt(x(1:end-1).^2 + x(2:end).^2))));
y(2)  = sum(abs(x).^(0.8) + 5 * sin(x.^3));
end

function F = kursaweFront()
% Obtained from jMetal
F = genetic.bench.loadFrontData('kursaweData')';
end