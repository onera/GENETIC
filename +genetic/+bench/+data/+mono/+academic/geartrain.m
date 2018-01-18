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
function bench = geartrain(n)
if nargin == 0
   bench = [];
   return
end
a        = 6.931;
% Output
bench    = struct('f'           ,@(x) geartrainFun(x,n,a),...
                  'bounds'      ,repmat([12 60],n,1),...
                  'dimsNumel'   ,repmat(length((12:60)),n,1));
%                                  'dsrchSpace'  ,mat2cell(repmat((12:60),n,1),ones(n,1), length((12:60))),...

end
% handle for computation
function y = geartrainFun(x,n,a)
y = (1/a - prod(x(1:floor(n/2)))/prod(x(floor(n/2)+1:end)))^2;
end
