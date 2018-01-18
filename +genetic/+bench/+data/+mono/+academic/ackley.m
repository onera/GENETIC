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
function bench = ackley(n)
if nargin == 0
   bench = [];
   return
end

a        = 20;
b        = 0.2;
c        = 2*pi;
% Output
bench    = struct('f'      , @(x) ackleyFun(x,n,a,b,c),...
                  'bounds' , repmat([-15 30],n,1),...
                  'xopt'   , zeros(n,1),...
                  'fopt'   , 0);

end
% handle for computation
function y = ackleyFun(x,n,a,b,c)
y = -a*exp(-b*sqrt(sum(x.^2)/n)) - exp(sum(cos(c*x))/n) + a + exp(1);
end
