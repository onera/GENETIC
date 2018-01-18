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
function bench = goldsteinprice()
n=2;
bench = struct('f'      , @(x) goldsteinpriceFun(x),...
               'bounds' , repmat([-2 2],n,1),...
               'fopt'   , 3,...
               'xopt'   , [0;-1],...
               'xDim'   , n);
end

function y = goldsteinpriceFun(x)
a =  (1 + (x(1) + x(2) + 1)^2*(19 - 14*x(1) + 3*x(1)^2 - 14*x(2) + 6*x(1)*x(2) + 3*x(2)^2));
b =  (30 + (2*x(1) - 3*x(2))^2*(18 - 32*x(1) + 12*x(1)^2 + 48*x(2) - 36*x(1)*x(2) + 27*x(2)^2));
y = a*b;
end
