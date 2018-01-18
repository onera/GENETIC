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
function bench = shaffer(n)
if nargin == 0
   bench = [];
   return
end
% SCHAFFER
a        = 0.5;
b        = 0.001;
bench    = struct('f'      , @(x) a + ( sin(sqrt(x'*x))^2 - a )/( (1 + b*(x'*x))^2 ),...
                  'bounds' , repmat([-100 100],n,1));
end
