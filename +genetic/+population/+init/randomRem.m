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
function x = randomRem(popSize, xMin, xMax, xDim)
% Output matrix <x> preallocation
x = zeros(xDim,popSize);
% Integer <nb> = # of random samples + + coordinates matrix <coordinatesMat>
n                    = 1e2*popSize;
estimatedMeanRadius  = 0.;
coordinatesMat       = rand(xDim,n);
% Estimated mean radius <estimatedMeanRadius> + distance <emd> computations
for j1 = 1:n-1
   for j2 = j1+1:n
      estimatedMeanRadius = estimatedMeanRadius + 2*sqrt((coordinatesMat(:,j1)-coordinatesMat(:,j2))'*(coordinatesMat(:,j1)-coordinatesMat(:,j2)))/(n*(n-1));
   end
end
% Square value of <estimatedMeanRadius> ---> <emd>
emd = estimatedMeanRadius*estimatedMeanRadius;
% Parameters initialisation integers <nok> + <nfo>
nok = 0; nfo = 0;
% MATLAB 'for' loop on <obj.popSize>
for j1 = 1:popSize
   flag = true;
   while flag
      x(:,j1) = rand(xDim,1);
      flag = false;
      for j2 = 1:j1-1
         d = (x(:,j1)-x(:,j2))'*(x(:,j1)-x(:,j2));
         if d < emd
            flag = true;
            nok = nok + 1;
            break
         end
      end
      if nok > j1*n
         nfo = nfo + 1;
         flag = false;
      end
   end
end
% Points dispatching within [xMin;xMax]
x = genetic.population.init.dispatchIndividuals(popSize, xMin, xMax, x);
end