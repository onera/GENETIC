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
function bench = cec17MaF12()
% MaF12 - CEC'2017 competition benchmark
%
%  This test problem is used to assess whether evolutionary many-objectives
%  algorithms are capable of dealing with scaled concave Pareto fronts with
%  complicated fitness landscapes.
%

% Bench parameters
M = 3;
K = M - 1;
L = 10;
%
n = K + L;
% if n ~= (K+L)
%    error('The benchmark ''cec17MaF12'' is only available for dimension %d.',K+L)
% end
% Output
bench    = struct('f'      , @(x) cec17MaF12Fun(x, M, K, L),...
                  'fDim'   , M,...
                  'bounds' , [zeros(n,1) (2:2:2*n)'],...
                  'fopt'   , [],...
                  'xDim'   , n);

end
% Handle for computation
function y = cec17MaF12Fun(x, M, K, L)
S              = (2:2:2*M)';
%
z              = x./(2:2:2*(K+L))';
%
t1             = zeros(K+L,1);
tmp            = flipud(cumsum(flipud(z(2:end,1)))./(1:K+L-1)');
expVector      = 0.02 + (50 - 0.02)*(0.98/49.98 - (1 - 2*tmp).*abs(floor(0.5 - tmp) + 0.98/49.98));
t1(1:K+L-1,1)  = z(1:K+L-1,1).^expVector;
t1(K+L,1)      = z(K+L,1);
%
t2             = zeros(K+L,1);
t2(1:K,1)      = 1 + (abs(t1(1:K,1) - 0.35) - 0.001).*(349.95*floor(t1(1:K,1) - 0.349)/0.349 + 649.95*floor(0.351 - t1(1:K,1))/0.649 + 1000);
t2(K+1:K+L,1)  = (1 + cos(122*pi*(0.5 - abs(t1(K+1:K+L,1) - 0.35)./(2*(floor(0.35 - t1(K+1:K+L,1)) + 0.35)))) + 380*(abs(t1(K+1:K+L,1) - 0.35)./(2*(floor(0.35 - t1(K+1:K+L,1)) + 0.35))).^2)/97;
%
t3             = zeros(M,1);
for i = 1:M-1
   is       = 1 + (i-1)*K/(M-1);
   ie       = i*K/(M-1);
   t3(i,1)  = r_nonsep(t2(is:ie),K/(M-1));
end
tmp = 0;
for i = K+1:K+L-1
   for j = i+1:K+L
      tmp = tmp + abs(t2(i,1) - t2(j,1));
   end
end
t3(M,1)        = (sum(t2(K+1:end,1),1) + tmp*2)/ceil(L/2)/(1 + 2*L - 2*ceil(L/2));
%
tmp            = zeros(M,1);
tmp(1:M-1,1)   = (t3(1:M-1,1) - 0.5)*max(1,t3(M,1)) + 0.5;
tmp(M,1)       = t3(M,1);
%
h              = flipud(cumprod([1;sin(tmp(1:end-1,1)*pi/2)],1)).*[1;cos(tmp(end-1:-1:1,1)*pi/2)];
%
y              = tmp(M,1)*ones(M,1) + S.*h;
end

% function F = cec17MaF12Front(M)
% F = genetic.bench.multi.loadData(['cec17MaF12Data' num2str(M) 'D'])';
% end

function vargout = r_nonsep(y, A)
    vargout = 0;
    for j = 1:size(y,1)
        tmp = 0;
        for k = 0:A-2
            tmp = tmp + abs(y(j,1) - y(1 + mod(j+k,size(y,1)),1));
        end
        vargout = vargout + y(j,1) + tmp;
    end
    vargout = vargout/(size(y,1)/A)/ceil(A/2)/(1 + 2*A - 2*ceil(A/2));
end