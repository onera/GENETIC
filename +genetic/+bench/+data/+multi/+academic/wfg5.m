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
function bench = wfg5(n)
if nargin == 0
   bench = [];
   return
end
% WFG5 benchmark
%
%  + ------------------------------------------------------ +
%  |  obj.  |  sep.  |  mod.  |     bia.     |     geo.     |
%  + ------------------------------------------------------ +
%  |  1..M  |   S    |   D    |       -      |   concave    |
%  + ------------------------------------------------------ +
%

% Bench parameters
M = 3;
K = M - 1;
%
if n < (K+1)
   warning('The benchmark ''wfg5'' dimension must be strictly greater or equal than %d.',K+1)
end
% Output
bench    = struct('f'      , @(x) wfg5Fun(x, M, K),...
                  'fDim'   , M,...
                  'bounds' , [zeros(n,1) (2:2:2*n)'],...
                  'fopt'   , wfg5Front(M));

end
% Handle for computation
function y = wfg5Fun(x, M, K)
n                 = length(x);
L                 = n - K;
S                 = (2:2:2*M)';
%
z                 = x./(2:2:2*n)';
%
t1                = shiftDeceptive(z, 0.35, 0.001, 0.05);
%
t2                = zeros(M,1);
weights           = ones(K/(M-1),1);
for i = 1:M-1
   is       = 1 + (i-1)*K/(M-1);
   ie       = i*K/(M-1);
   t2(i,1)  = weightedSum(t1(is:ie,1), weights);
end
is                = K + 1;
ie                = K + L;
weights           = ones(L,1);
t2(M,1)           = weightedSum(t1(is:ie,1), weights);
%
tmp               = zeros(M,1);
tmp(1:M-1,1)      = (t2(1:M-1,1) - 0.5)*max(1,t2(M,1)) + 0.5;
tmp(M,1)          = t2(M,1);
%
h                 = flipud(cumprod([1;sin(tmp(1:M-1,1)*pi/2)],1)).*[1;cos(tmp(M-1:-1:1,1)*pi/2)];
%
y                 = tmp(M,1)*ones(M,1) + S.*h;
end

function F = wfg5Front(M)
F = genetic.bench.loadFrontData(['wfg5Data' num2str(M) 'D'])';
end

function vargout = shiftDeceptive(y, A, B, C)
lambda   = (1 - C + (A-B)/B)/(A-B);
mu       = (1 - C + (1-A-B)/B)/(1-A-B);
vargout  = 1 + (abs(y - A) - B).*(lambda*floor(y - (A-B)) + mu*floor((A+B) - y) + 1/B);
end

function vargout = weightedSum(y, w)
vargout = sum(w.*y,1)/sum(w,1);
end