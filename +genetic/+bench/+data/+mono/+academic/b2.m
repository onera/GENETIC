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
function bench = b2()
n        = 2;
a        = 2;
b        = 0.3;
c        = 0.4;
d        = 0.7;
bench    = struct('f'      , @(x) x(1)^2 + a*x(2)^2 - b*cos(3*pi*x(1)) - c*cos(4*pi*x(2)) + d,...
                  'bounds' , repmat([-2 2],n,1));
%
bench.xDim = n;
end
