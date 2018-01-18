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
function bench = cec17MaF11()
% MaF11 - CEC'2017 competition benchmark
%
%  This test problem is used to assess whether evolutionary many-objectives
%  algorithms are capable of dealing with scaled disconnected Pareto fronts.
%

% Bench parameters
M = 3;
K = M - 1;
L = 10;
%
n = K + L;
% if n ~= (K+L)
%    error('The benchmark ''cec17MaF11'' is only available for dimension %d.',K+L)
% end
% Output
bench    = struct('f'      , @(x) cec17MaF11Fun(x, M, K, L),...
                  'fDim'   , M,...
                  'bounds' , [zeros(n,1) (2:2:2*n)'],...
                  'fopt'   , [],...
                  'xDim'   , n);

end
% Handle for computation
function y = cec17MaF11Fun(x, M, K, L)
S                 = (2:2:2*M)';
%
z                 = x./(2:2:2*(K+L))';
%
t1                = zeros(K+L,1);
t1(1:K,1)         = z(1:K,1);
t1(K+1:end,1)     = abs(z(K+1:end,1) - 0.35)./(abs(floor(0.35 - z(K+1:end,1))) + 0.35);
%
t2                = zeros(K+L,1);
t2(1:K,1)         = t1(1:K,1);
t2(K+1:K+L/2,1)   = t1(K+1:2:end-1,1) + t1(K+2:2:end,1) + 2*abs(t1(K+1:2:end-1,1) - t1(K+2:2:end,1));
t2(K+1:end,1)     = genetic.tools.roundn(t2(K+1:end,1),-6);
%
t3                = zeros(M,1);
for i = 1:M-1
   is       = 1 + (i-1)*K/(M-1);
   ie       = i*K/(M-1);
   weights  = ones(K/(M-1),1);
   t3(i,1)  = weightedSum(t2(is:ie,1), weights);
end
is                = K + 1;
ie                = K + L/2;
weights           = ones(L/2,1);
t3(M,1)           = weightedSum(t2(is:ie,1), weights);
%
tmp               = zeros(M,1);
tmp(1:M-1,1)      = (t3(1:M-1,1) - 0.5)*max(1,t3(M,1)) + 0.5;
tmp(M,1)          = t3(M,1);
%
h                 = flipud(cumprod([1;1-cos(tmp(1:M-1,1)*pi/2)],1)).*[1;1-sin(tmp(M-1:-1:1,1)*pi/2)];
h(M,1)            = 1 - tmp(1,1)*(cos(5*pi*tmp(1,1)))^2;
%
y                 = tmp(M,1)*ones(M,1) + S.*h;
end

% function F = cec17MaF11Front(M)
% F = genetic.bench.multi.loadData(['cec17MaF11Data' num2str(M) 'D'])';
% end

function vargout = weightedSum(y, w)
vargout = sum(w.*y,1)/sum(w,1);
end