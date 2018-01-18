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
classdef simpleScheme < handle
   properties
      nGen = 0;
   end
   methods
      function self = simpleScheme()
      end
      %%
      function [xopt, fopt, info] = apply(self, x0)
         x0Provided = true;
         if nargin < 2 || isempty(x0)
            x0          = [];
            x0Provided  = false;
         end
         self.printHeader(x0Provided);
         self.printStart();
         % Initialisation -------------------------------------------------
         % Group creation
         group = self.makeGroup();
         % Initial population
         X0    = self.initPopulation(x0);
         group.moveTo(X0);
         % Initial evaluation
         Y0    = self.simulator.eval(group);
         group.tellObjective(Y0)
         self.postEval(group);
         %
         self.printIterHead()
         self.printIter(group);
         % Main optimization loop -----------------------------------------
         while ~self.stop(group)
            self.nGen   = self.nGen + 1;
            group       = self.evolve(group);
            Y           = self.simulator.eval(group);
            group.tellObjective(Y);
            self.postEval(group);
            %
            if mod(self.printCtr, self.modPrintHead) == 0
               self.printCtr = 1;
               self.printIterHead();
            end
            self.printIter(group);
         end
         self.printEnd();
         % Output variable assignement ------------------------------------
         [xopt, fopt] = group.getAbsoluteBest();
         %
         info = self.fillInfo(group);
      end
      %% 
      function postEval(self, group)
      end
   end
end
