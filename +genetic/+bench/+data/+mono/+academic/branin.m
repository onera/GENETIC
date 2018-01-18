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
function bench = branin()
n = 2;
% BRANIN
a        = 1;
b        = 5.1/(4*pi^2);
c        = 5/pi;
d        = 6;
e        = 10;
f        = 1/(8*pi);
bench    = struct('f'      , @(x) a*(x(2) - b*x(1)^2 + c*x(1) - d)^2 + e*(1 - f)*cos(x(1)) + e,...
                  'bounds' , [-5 10;0 15],...
                  'fopt'   , 0.397887,...
                  'xopt'   , [[-pi;12.275],[pi;2.275],[9.42478;2.475]],...
                  'xDim'   , n);

end
