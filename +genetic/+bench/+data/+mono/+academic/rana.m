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
function bench = rana(n)
if nargin == 0
   bench = [];
   return
end
% RANA
bench = struct('f'      , @(x) ranaFun(x, n),...
               'bounds' , repmat([-520 520],n,1));
end

function y = ranaFun(x, n)
y = 0;
for i=1:n
   if i==1
      xr = x(2);
   elseif i==2,
      xr = x(1);
   end
   y = y + x(i)*sin(sqrt(abs(xr+1-x(i))))*cos(sqrt(abs(xr+1+x(i)))) + (xr+1)*cos(sqrt(abs(xr+1-x(i))))*sin(sqrt(abs(xr+1+x(i))));
end
y = y/n;
end
