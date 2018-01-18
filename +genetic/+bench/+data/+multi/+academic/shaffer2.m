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
function bench = shaffer2()
n = 1;
% Output
bench    = struct('f'      , @(x) shaffer2Fun(x),...
                  'fDim'   , 2,...
                  'bounds' , repmat([-5 10],n,1),...
                  'fopt'   , [],...
                  'xDim'   , n);

end
% handle for computation
function y = shaffer2Fun(x)
y     = zeros(2,1);
if x <= 1
   y(1) = -x;
elseif x <= 3 && x > 1
   y(1)  = x-2;
elseif 3 < x && x <= 4
   y(1) = 4-x;
else
   y(1) = x-4;
end
y(2)  = (x-5)^2;
end
