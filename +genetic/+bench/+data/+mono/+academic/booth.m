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
function bench = booth()
n=2;
%
a        = 2;
b        = 7;
c        = 5;
bench    = struct('f'      , @(x) (x(1) + a*x(2) - b)^2 + (a*x(1) + x(2) - c)^2,...
                  'bounds' , repmat([-10 10],n,1),...
                  'xopt'   , [1;3],...
                  'fopt'   , 0,...
                  'xDim'   , n);
end
