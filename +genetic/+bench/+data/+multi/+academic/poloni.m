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
function bench = poloni()
n = 2;
% Output
bench    = struct('f'      , @(x) poloniFun(x),...
                  'fDim'   , 2,...
                  'bounds' , repmat([-pi pi],n,1),...
                  'fopt'   , [],...
                  'xDim'   , n);

end
% handle for computation
function y = poloniFun(x)
y     = zeros(2,1);
A1    = 0.5 * sin(1) - 2 * cos(1) + sin(2) - 1.5 * cos(2);
A2    = 1.5 * sin(1) - cos(1) + 2 * sin(2) - 0.5 * cos(2);
B1    = 0.5 * sin(x(1)) - 2 * cos(x(1)) + sin(x(2)) - 1.5 * cos(x(2));
B2    = 1.5 * sin(x(1)) - cos(x(1)) + 2 * sin(x(2)) - 0.5 * cos(x(2));
y(1)  = 1 + (A1 - B1)^2 + (A2 - B2)^2;
y(2)  = (x(1) + 3)^2 + (x(2) + 1)^2;
end

% function F = poloniFront()
% % Obtained from jMetal
% F = genetic.bench.multi.loadData('poloniData')';
% end