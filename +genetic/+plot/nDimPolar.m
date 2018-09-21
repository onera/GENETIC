% See the article:
%   W.L. Goffe, Visualizing Multi-Dimensional Functions in
%   Economics. 1999
        
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
function [t,X] = nDimPolar(n, N)
% Real parameter for the angle
t               = linspace(0,2 * pi, N);

% Generalised Euler angles
theta           = zeros(n-1,N);
theta(n-1,:)    = t;
for i =1:n-2
    theta(i,:)  = pi/2 * sin(2^(n - i -1) * t);
end

% Coordinates
X = ones(n,N);
for i = 1:n
    for j = 1:n-i+1
        if j < n-i+1
            X(i,:) = X(i,:) .* cos(theta(j,:));
        else
            if j <= n-1
                X(i,:) = X(i,:) .* sin(theta(j,:));
            end
        end
    end
end
% 

end