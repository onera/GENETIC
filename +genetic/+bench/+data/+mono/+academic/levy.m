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
function bench = levy(n)
if nargin == 0
   bench = [];
   return
end
% LEVY
bench = struct('f'      , @(x)levyFun(x,n),...
               'bounds' , repmat([-10 10],n,1),...
               'xopt'   , ones(n,1),...
               'fopt'   , 0);
end

function y = levyFun(x,n)
x = 1 + (x-1)/4;
y = (sin(pi*x(1)))^2 + sum(((x(1:n-1) - 1).^2).*(1 + 10*(sin(pi*x(1:n-1) + 1)).^2)) + (x(n) - 1)^2*(1 + 10*sin(2*pi*x(n))^2);
end
