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
function bench = cec17MaF14()
% MaF14 - CEC'2017 competition benchmark
%
%  This test problem is used to assess whether evolutionary many-objectives
%  algorithms are capable of dealing with complicated fitness landscape with
%  mixed variable separability, especially in large-scale cases.
%

% Bench parameters
M        = 3;
Nk       = 2;
c        = 3.8*0.1*(1 - 0.1);
%
for i = 1:M-1
   c = [c;3.8*c(end,1)*(1 - c(end,1))];
end
%
S        = ceil(round(c./sum(c)*(20*M))/Nk);
L        = [0;cumsum(Nk*S)];
n = M - 1 + L(end);
% if n ~= (M-1)+L(end)
%    error('The benchmark ''cec17MaF14'' is only available for dimension equal to %d.',(M-1)+L(end))
% end
% Output
bench    = struct('f'      , @(x) cec17MaF14Fun(x, M, Nk, S, L),...
                  'fDim'   , M,...
                  'bounds' , [repmat([0 1],M-1,1);repmat([0 10],L(end),1)],...
                  'fopt'   , [],...
                  'xDim'   , n);

end
% Handle for computation
function y = cec17MaF14Fun(x, M, Nk, S, L)
n           = M - 1 + L(end);
tmp         = x;
tmp(M:n,1)  = (1 + (M:n)'/n).*tmp(M:n,1) - 10*tmp(1,1)*ones(n-M+1,1);
g           = zeros(M,1);
%
for i = 1:2:M
   for k = 1:Nk
      is       = L(i) + M - 1 + (k - 1)*S(i) + 1;
      ie       = L(i) + M - 1 + k*S(i);
      g(i,1)   = g(i,1) + fun1Local(tmp(is:ie,1));
   end
end
%
for i = 2:2:M
   for k = 1:Nk
      is       = L(i) + M - 1 + (k - 1)*S(i) + 1;
      ie       = L(i) + M - 1 +k*S(i);
      g(i,1)   = g(i,1) + fun2Local(tmp(is:ie,1));
   end
end
%
g           = g./S/Nk;
y           = (1 + g).*flipud(cumprod([1;tmp(1:M-1,1)])).*[1;1-tmp(M-1:-1:1,1)];
end

% function F = cec17MaF14Front(M)
% F = genetic.bench.multi.loadData(['cec17MaF14Data' num2str(M) 'D'])';
% end

function y = fun1Local(x)
y = sum(x.^2 - 10*cos(2*pi*x) + 10,1);
end

function y = fun2Local(x)
y = sum(100*(x(1:length(x)-1,1).^2 - x(2:length(x),1)).^2 + (x(1:length(x)-1,1) - 1).^2,1);
end