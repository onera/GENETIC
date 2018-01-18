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
function bench = bohachevsky2()
n = 2;
% BOHACHEVSKY-2
a        = 1;
b        = 2;
c        = 0.3;
w1       = 3*pi;
w2       = 4*pi;
bench    = struct('f'      , @(x) a*x(1)^2 + b*x(2)^2 - c*cos(w1*x(1))*cos(w2*x(2)) + c,...
                  'bounds' , repmat([-100 100],n,1),...
                  'xopt'   , [0;0],...
                  'fopt'   , 0,...
                  'xDim'   , n);
             
end
