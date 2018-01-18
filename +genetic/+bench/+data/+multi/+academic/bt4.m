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
function bench = bt4()
% BT4 benchmark - Multi-objective test problems with bias
%
%  This multi-objective test problem  is characterized by two criteria that
%  must be optimized. The number of decision variables is set to 30.
%

% Bench parameters
n = 30;
M = 2;

% Output
bench    = struct('f'      , @(x) bt4Fun(x, M),...
                  'fDim'   , M,...
                  'bounds' , repmat([0 1],n,1),...
                  'fopt'   , [],...
                  'xDim'   , n);

end
% Handle for computation
function y = bt4Fun(x, M)

I1         = 2:2:D;
I2         = 3:2:D;
Y          = X - sin(repmat(1:D,N,1)*pi/2/D);
temp1      = X(:,1) < 0.25;
temp2      = 0.25 <= X(:,1) & X(:,1) < 0.5;
temp3      = 0.5 <= X(:,1) & X(:,1) < 0.75;
temp4      = 0.75 <= X(:,1);
X(temp1,1) = (1-(1-4*X(temp1,1)).^0.06)/4;
X(temp2,1) = (1+(4*X(temp2,1)-1).^0.06)/4;
X(temp3,1) = (3-(3-4*X(temp3,1)).^0.06)/4;
X(temp4,1) = (3+(4*X(temp4,1)-3).^0.06)/4;
PopObj(:,1) = X(:,1)         + sum(Y(:,I1).^2+(1-exp(-Y(:,I1).^2/1e-8))/5,2);
PopObj(:,2) = 1-sqrt(X(:,1)) + sum(Y(:,I2).^2+(1-exp(-Y(:,I2).^2/1e-8))/5,2);

y     = zeros(M,1);
n     = length(x);
ind1  = (2:2:n);
ind2  = (3:2:n);
tmp   = x - sin(pi*(1:n)'/2/n);
tmp1  = x(1,1) < 0.25;
tmp2  = 0.25 <= x(1,1) & x(1,1) < 0.5;
tmp3  = 0.5 <= x(1,1) & x(1,1) < 0.75;
tmp4  = 0.75 <= x(1,1);

%% !! ICI !! %%

y(1)  = abs(x(1,1))^0.02 + sum(tmp(ind1,1).^2 + (1 - exp(-tmp(ind1,1).^2/1e-8))/5,1);
y(2)  = 1 - sqrt(abs(x(1,1))^0.02) + sum(tmp(ind2,1).^2 + (1 - exp(-tmp(ind2,1).^2/1e-8))/5,1);
end

% function F = bt4Front(M)
% F = genetic.bench.multi.loadData(['bt4Data' num2str(M) 'D'])';
% end