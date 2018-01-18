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
function bench = dtlz2()
M = 3;
K = 10;
n = M-1 + K;

% Output
bench    = struct('f'      , @(x) dtlz2Fun(x, M, K),...
                  'fDim'   , M,...
                  'bounds' , repmat([0 1],n,1),...
                  'fopt'   , dtlz2Front(M),...
                  'xDim'   , n);

end
% handle for computation
function y = dtlz2Fun(x, M, K)
n           = (M - 1) + K;
xm          = x(n-K+1:end,1);
g           = sum((xm - 0.5).^2,1);
y(1,1)      = (1 + g).*prod(cos(pi/2*x(1:M-1,:)),1);
for i = 2:M-1
   y(i,1) = (1 + g).*prod(cos(pi/2*x(1:M-i,:)),1).*sin(pi/2*x(M-i+1,:));
end
y(M,1)      = (1 + g).*sin(pi/2*x(1,1));
end

%         xm = x(n-k+1:end,:); %xm contains the last k variables
%         g = sum((xm - 0.5).^2, 1);
% 
%         % Computes the functions
%         f(1,:) = (1 + g).*prod(cos(pi/2*x(1:M-1,:)),1);
%         for ii = 2:M-1
%            f(ii,:) = (1 + g) .* prod(cos(pi/2*x(1:M-ii,:)),1) .* ...
%               sin(pi/2*x(M-ii+1,:));
%         end
%         f(M,:) = (1 + g).*sin(pi/2*x(1,:));

function F = dtlz2Front(M)
% Obtained from jMetal
F = genetic.bench.loadFrontData(['dtlz2Data' num2str(M) 'D'])';
end