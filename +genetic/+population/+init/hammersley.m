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
function x = hammersley(sizePop, xMin, xMax, xDim)
% Output matrix <x> preallocation
x = zeros(xDim,sizePop);
% 
MaxPrime = 1600;
job      = genetic.population.init.prime(1:MaxPrime);
% Dimension # random selection ---> integer <dimSelect>
dimSelect = genetic.tools.randomSamples((1:xDim));
%
for j = 1:sizePop
    x(dimSelect,j) = xMin(dimSelect) + abs(xMax(dimSelect)-xMin(dimSelect))*(j-1)/(sizePop-1);
    k = dimSelect;
    % 
    for i = 1:xDim-1
        k = k + 1;
        if k > xDim
            k = 1;
        end
        nbp   = job(k);
        pstar = nbp;
        kstar = j;
        phi   = 0;
        % 
        while kstar > 0
            a     = mod(kstar,nbp);
            phi   = phi + a/pstar;
            kstar = fix(kstar/nbp);
            pstar = pstar*nbp;
        end
        x(k,j) = xMin(k) + phi*abs(xMax(k)-xMin(k));
    end
end
        
end
