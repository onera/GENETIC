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
function bench = wfg3(n)
if nargin == 0
   bench = [];
   return
end
% WFG3 benchmark
%
%  + ------------------------------------------------------ +
%  |  obj.  |  sep.  |  mod.  |     bia.     |     geo.     |
%  + ------------------------------------------------------ +
%  |  1..M  |   NS   |   U    |       -      |   linear,    |
%  |        |        |        |              | degenerate   |
%  + ------------------------------------------------------ +
%

% Bench parameters
M = 3;
K = M - 1;
%
if n < (K+1)
   warning('The benchmark ''wfg3'' dimension must be strictly greater or equal than %d.',K+1)
end
% Output
bench    = struct('f'      , @(x) wfg3Fun(x, M, K),...
                  'fDim'   , M,...
                  'bounds' , [zeros(n,1) (2:2:2*n)'],...
                  'fopt'   , wfg3Front(M));

end
% Handle for computation
function y = wfg3Fun(x, M, K)
n                 = length(x);
L                 = n - K;
S                 = (2:2:2*M)';
%
z                 = x./(2:2:2*n)';
%
t1                = zeros(n,1);
t1(1:K,1)         = z(1:K,1);
t1(K+1:end,1)     = abs(z(K+1:end,1) - 0.35)./abs(floor(0.35 - z(K+1:end,1)) + 0.35);
%
t2                = zeros(K+L/2,1);
t2(1:K,1)         = t1(1:K,1);
t2(K+1:K+L/2,1)   = (t1(K+1:2:end,1) + t1(K+2:2:end,1) + 2*abs(t1(K+1:2:end,1) - t1(K+2:2:end,1)))/3;
%
t3                = zeros(M,1);
weights           = ones(K/(M-1),1);
for i = 1:M-1
   is       = 1 + (i-1)*K/(M-1);
   ie       = i*K/(M-1);
   t3(i,1)  = weightedSum(t2(is:ie,1), weights);
end
is                = K + 1;
ie                = K + L/2;
weights           = ones(L/2,1);
t3(M,1)           = weightedSum(t2(is:ie,1), weights);
%
tmp               = zeros(M,1);
tmp(1:M-1,1)      = (t3(1:M-1,1) - 0.5).*max([[1;zeros(M-2,1)] t3(M,1)*ones(M-1,1)],[],2) + 0.5;
tmp(M,1)          = t3(M,1);
%
h                 = flipud(cumprod([1;tmp(1:M-1,1)],1)).*[1;1-tmp(M-1:-1:1,1)];
%
y                 = tmp(M,1)*ones(M,1) + S.*h;
end

function F = wfg3Front(M)
F = genetic.bench.loadFrontData(['wfg3Data' num2str(M) 'D'])';
end

function vargout = weightedSum(y, w)
vargout = sum(w.*y,1)/sum(w,1);
end