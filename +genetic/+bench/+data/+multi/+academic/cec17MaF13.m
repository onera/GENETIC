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
function bench = cec17MaF13(n)
if nargin == 0
   bench = [];
   return
end
% MaF13 - CEC'2017 competition benchmark
%
%  This test problem is used to assess whether evolutionary many-objectives
%  algorithms are capable of dealing with degenerate Pareto fronts and complex
%  variable linkages.
%

% Bench parameters
M = 3;
if n < 2
   error('The benchmark ''cec17MaF13'' is only available for dimension greater or equal to 2.')
end
% Output
bench    = struct('f'      , @(x) cec17MaF13Fun(x, M),...
                  'fDim'   , M,...
                  'bounds' , [repmat([0 1],2,1);repmat([-2 2],n-2,1)],...
                  'fopt'   , []);

end
% Handle for computation
function y = cec17MaF13Fun(x, M)
n        = length(x);
tmp      = x - 2*x(2,1)*ones(n,1).*sin(2*pi*x(1,1)*ones(n,1) + pi*(1:n)'/n);
y        = zeros(M,1);
y(1,1)   = sin(x(1,1)*pi/2) + 2*mean(tmp(4:3:n,1).^2);
y(2,1)   = cos(x(1,1)*pi/2)*sin(x(2,1)*pi/2) + 2*mean(tmp(5:3:n,1).^2);
y(3,1)   = cos(x(1,1)*pi/2)*cos(x(2,1)*pi/2) + 2*mean(tmp(3:3:n,1).^2);
y(4:M,1) = repmat(y(1,1)^2 + y(2,1)^10 + y(3,1)^10 + 2*mean(tmp(4:n,1).^2),M-3,1);
end

% function F = cec17MaF13Front(M)
% F = genetic.bench.multi.loadData(['cec17MaF13Data' num2str(M) 'D'])';
% end