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
function bench = ellipsoid(n)
if nargin == 0
   bench = [];
   return
end
% ELLIPSOID
% x(1)^2 + (x(1) + x(2))^2 + ... + (x(1) + ... + x(n))^2
M        = triu(ones(n,n));
bench    = struct('f'      , @(x) sum((x'*M).^2),...
                  'bounds' , repmat([-65.536 65.536],n,1));
end
