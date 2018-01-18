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
function x = regularise(x, sizePop, xMin, xMax, xDim)


if sizePop == 1
    return
else
    for i = 1:xDim
        if i == 1
            job = genetic.population.init.hammersley(sizePop, xMin(i), xMax(i), 1);
            x(i,:) = job(randperm(sizePop));
        else
            equalFlag = false;
            for j = i-1:-1:1
                if xMin(i) == xMin(j) && xMax(i) == xMax(j)
                    x(i,:) = x(j,randperm(sizePop));
                    equalFlag = true;
                    break
                end
            end
            if ~equalFlag
                job = genetic.population.init.hammersley(sizePop, xMin(i), xMax(i), 1);
                x(i,:) = job(randperm(sizePop));
            end
        end
    end
end

end