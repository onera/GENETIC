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
classdef multi < genetic.optimizer.base
   properties
       metric = 'hvmc';
       tolY     = 1e-8;
      nTolY_   = 0;
      nTolY    = 1;

   end
   methods
      function self = multi(varargin)
         self@genetic.optimizer.base(varargin{:});
         % specific stopping criteria
         if ~isempty(self.metric)
            yTolReached    = struct('name', 'yTolReached', 'cleanStop', true, 'condByConstraints', true);
         end
         self.addStopTest(yTolReached);

         % Things to display
         imp            = struct('name','I','call','printImp','align','c','dim',1);
         nGen           = struct('name','nGen','call','printnGen','align','r','dim',6);
         neval          = struct('name','nEval/maxFunEval','call','printNeval','align','r','dim',16);
         self.addDisplayElement(imp, nGen, neval);

      end
      %% yTolReached
      function [out, msg] = yTolReached(obj, group)
         out   = false;
         msg   = '';         
         if group.hasImproved
            if ~isempty(group.bestObjectiveChange)
               out = group.bestObjectiveChange < obj.tolY;
            end
            if out
               obj.nTolY_     = obj.nTolY_ + 1;
               if obj.nTolY_ == obj.nTolY
                  msg            = sprintf('Best objective decrease is below tolerance: %1.2e (tol: %1.2e)',group.bestObjectiveChange, obj.tolY);
               else
                  out = false;
               end
            else
               obj.nTolY_ = 0;
            end
         end
      end
            % peut etre a deplacer
      function group = makeGroup(self)
         group = genetic.population.mgroup(self.popSize);
         try
            assignTopo = true;
            group.initTopology(self.topologyType, assignTopo);
         catch
         end
      end
   end
   %%
   methods (Static)
      function mx = getFixedRowSumIntegerMatrix(nObj, rowSum)
         if nObj == 1
            mx = rowSum; 
            return;
         end
         %
         mx = [];
         for k = 0:rowSum
            tmp   = genetic.optimizer.multi.getFixedRowSumIntegerMatrix(nObj-1, rowSum-k);
            mx    = [mx;k*ones(size(tmp,1),1) tmp];
         end
         %
      end
   end
end