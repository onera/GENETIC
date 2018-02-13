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
classdef mono < genetic.optimizer.base
   properties
      tolX     = 1e-8;
      tolY     = 1e-8;
      targetY  = -inf;
      %
      nTolX    = 1;
      nTolY    = 1;
      nTolX_   = 0;
      nTolY_   = 0;
      %
   end
   methods
      function self = mono(varargin)
         self@genetic.optimizer.base(varargin{:});
         % Adding new stopping criteria
         yTolReached    = struct('name', 'yTolReached', 'cleanStop', true, 'condByConstraints', true);
         xTolReached    = struct('name', 'xTolReached', 'cleanStop', true, 'condByConstraints', true);
         yTargetReached = struct('name', 'yTargetReached', 'cleanStop', true, 'condByConstraints', true);
         self.addStopTest(yTolReached, xTolReached, yTargetReached);
         % things to print during the iteration
         imp            = struct('name','I','call','printImp','align','c','dim',1);
         nGen           = struct('name','nGen','call','printnGen','align','r','dim',6);
         bcf            = struct('name','fGen*','call','printBestCurrentF','align','c','dim',16);
         bf             = struct('name','f*','call','printBestF','align','c','dim',16);
         neval          = struct('name','nEval/maxFunEval','call','printNeval','align','r','dim',16);
         ytol           = struct('name','funDecr', 'call','printTolY','align','c','dim',16);
         xtol           = struct('name','xMove', 'call','printTolX','align','c','dim',16);
         self.addDisplayElement(imp, nGen, bcf, bf, neval, ytol, xtol);
         if ~self.isUnconstrained
            cstV           = struct('name','cstViolation','call','printCst','align','c','dim',16);
            self.addDisplayElement(cstV);
         end
         %
         if self.printOnlyImp
            self.removeDisplay('fGen*');
         end
         %
      end
      %
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
      %% xTolReached
      function [out,msg] = xTolReached(obj, group)
         out   = false;
         msg   = '';
         if group.hasImproved
            if ~isempty(group.bestObjectiveChange)
               out = group.bestValueChange < obj.tolX;
            end
            if out
               obj.nTolX_     = obj.nTolX_ +1;
               if obj.nTolX_ == obj.nTolX
                  msg = sprintf('Best point change is below tolerance: %1.2e (tol: %1.2e)',group.bestValueChange, obj.tolX);
               else
                  out = false;
               end
            else
               obj.nTolX_ = 0;
            end
         end
      end
      %% yTargetReached
      function [out, msg] = yTargetReached(obj, group)
         out = false;
         if ~isempty(group.best)
            out = group.best.objective <= obj.targetY;
         end
         msg = '';
         if out
            msg = sprintf('Best objective below target: %1.2e (target: %1.2e)',group.best.objective,obj.targetY);
         end
      end

      %%
      % peut etre a deplacer
      function group = makeGroup(self)
         group = genetic.population.group(self.popSize);
         try
            assignTopo = true;
            group.initTopology(self.topologyType, assignTopo);
         catch
         end
      end
      %%
      function out = printBestF(self, group)
         out = sprintf('%1.5e',group.best.objective);
      end
      function out = printBestCurrentF(self, group)
         out = sprintf('%1.5e', min(group.getObjective()));
      end

      function out = printFlag(self, group)
         out            = self.iterFlag;
         self.iterFlag  = '';
      end

      function out = printTolX(self, group)
         if group.hasImproved && self.nGen > 0
            out = sprintf('%1.5e',group.bestValueChange);
         else
            out = sprintf('-');
         end
      end
      function out = printTolY(self, group)
         if group.hasImproved && self.nGen > 0
            out = sprintf('%1.5e', group.bestObjectiveChange);
         else
            out = sprintf('-');
         end
      end
      %%
      function out = printCst(self, group)
         if group.hasImproved
            out = sprintf('%1.5e', self.constraints.maxViolation(group.best.value, false));
         else
            out = sprintf('-');
         end
      end
      %
      %%
   end
end