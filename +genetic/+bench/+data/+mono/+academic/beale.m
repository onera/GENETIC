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
function bench = beale()

n        = 2;
a        = 1.5;
b        = 2.25;
c        = 2.625;
bench    = struct('f'      ,@(x) (a - x(1)*(1 - x(2)))^2 + (b - x(1)*(1 - x(2)^2))^2 + (c - x(1)*(1 - x(2)^3))^2,...
                  'bounds' ,repmat([-4.5 4.5],n,1),...
                  'fopt'   , 0,...
                  'xopt'   , [3;0.5],...
                  'xDim'   , n );
end
