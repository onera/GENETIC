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
function bench = rastrigin(n)
if nargin == 0
   bench = [];
   return
end
% RASTRIGIN
a        = 10;
b        = 2*pi;
bench    = struct('f'      , @(x) a*n + sum(x.^2 - a*cos(b*x)),...
                  'bounds' , repmat([-5.12 5.12],n,1));
end
