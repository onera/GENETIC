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
function bench = wfg6(n)
if nargin == 0
   bench = [];
   return
end
% WFG6 benchmark
%
%  + ------------------------------------------------------ +
%  |  obj.  |  sep.  |  mod.  |     bia.     |     geo.     |
%  + ------------------------------------------------------ +
%  |  1..M  |   NS   |   U    |       -      |   concave    |
%  + ------------------------------------------------------ +
%

% Bench parameters
M = 3;
K = M - 1;
%
if n < (K+1)
   warning('The benchmark ''wfg6'' dimension must be strictly greater or equal than %d.',K+1)
end
% Output
bench    = struct('f'      , @(x) wfg6Fun(x, M, K),...
                  'fDim'   , M,...
                  'bounds' , [zeros(n,1) (2:2:2*n)'],...
                  'fopt'   , wfg6Front(M));

end
% Handle for computation
function y = wfg6Fun(x, M, K)
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
t2                = zeros(M,1);
for i = 1:M-1
   is       = 1 + (i-1)*K/(M-1);
   ie       = i*K/(M-1);
   t2(i,1)  = nonSepReduction(t1(is:ie,1), K/(M-1));
end
job               = 0;
for i1 = K+1:K+L
   for i2 = i1+1:K+L
      job = job + abs(t1(i1,1) - t1(i2,1));
   end
end
t2(M,1)           = (sum(t1(K+1:end,1),1) + 2*job)/ceil(L/2)/(1 + 2*L - 2*ceil(L/2));
%
tmp               = zeros(M,1);
tmp(1:M-1,1)      = (t2(1:M-1,1) - 0.5)*max(1,t2(M,1)) + 0.5;
tmp(M,1)          = t2(M,1);
%
h                 = flipud(cumprod([1;sin(tmp(1:M-1,1)*pi/2)],1)).*[1;cos(tmp(M-1:-1:1,1)*pi/2)];
%
y                 = tmp(M,1)*ones(M,1) + S.*h;
end

function F = wfg6Front(M)
F = genetic.bench.loadFrontData(['wfg6Data' num2str(M) 'D'])';
end

function vargout = nonSepReduction(y, A)
vargout = 0;
for i = 1:length(y)
   tmp = 0;
   for k = 0:A-2
      tmp = tmp + abs(y(i,1) - y(1 + mod(i+k,length(y)),1));
   end
   vargout = vargout + y(i,1) + tmp;
end
vargout = vargout/(length(y)/A)/ceil(A/2)/(1+2*A-2*ceil(A/2));
end