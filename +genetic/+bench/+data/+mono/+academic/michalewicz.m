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
function bench = michalewicz(n)
if nargin == 0
   bench = [];
   return
end
% MICHALEWICZ
m        = 10;
bench    = struct('f'      , @(x) michalewiczFun(x,m),...
                  'bounds' , repmat([0 pi],n,1));
end

function y = michalewiczFun(x,m)
s = 0;
for i=1:n
   s = s + sin(x(i))*(sin(i*x(i)^2/pi))^(2*m); 
end; 
y = 9.66 - s;
end
