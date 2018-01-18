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
function [benchExists, path] = checkAvailability(benchName)
% Checks whether a benchmark exists based on its name.
[pathToBench,EXT] = genetic.bench.getPath();
benchExists       = false;
for i = 1:length(pathToBench)
   p        = pathToBench{i};
   if exist([p,benchName,EXT{i}],'file') > 0
      benchExists = true;
      path        = [p,benchName];
      break
   end
end
end