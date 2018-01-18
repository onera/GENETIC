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
classdef mindividual < genetic.population.individual
   methods
      %% Constructor
      function obj = mindividual()
         obj@genetic.population.individual()
      end      
      %% tellObjective
      function tellObjective(obj, newObjective)
%          obj.isWaitingObjective  = false;
         obj.previousObjective   = obj.objective;
         obj.objective           = newObjective;
         obj.best.value          = obj.value;
         obj.best.objective      = newObjective;
         % Submit the new objective to the mother group to test if it
         % belongs to the pareto front
         obj.mother.updateMemory(obj.best)
         %
         if ~obj.isAlive
            obj.kill();
         end
      end
      
   end
end