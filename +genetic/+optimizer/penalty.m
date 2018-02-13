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
classdef penalty < handle%genetic.optimizer.base
   properties
      nPen           = 0;
      maxPen         = inf;
      sigma          = 1;
      sigmaMax       = 1e8;
   end
   properties (SetAccess = protected)
      initFun        = [];
      initOpt        = [];
      initCst        = [];
      constraints    = [];
      toPenalise     = {};
      subOptimizer   = [];
   end
   %
   methods
      function self = penalty(subOptimizer, toPenalise, constraints, options)
         self.constraints  = genetic.constraints(subOptimizer.constraints);
         self.initOpt      = options;
         self.initCst      = constraints;
         self.initFun      = subOptimizer.simulator.f_;
         %
         self.subOptimizer = subOptimizer;
         self.toPenalise   = toPenalise;
      end
      %%
      function [xopt, fopt, info] = apply(self, x0)
         self.subOptimizer.constraints.skip  = self.toPenalise;
         [penF, d]                           = self.subOptimizer.constraints.penalise(self.toPenalise);
         % Initial unconstrained subproblem
%          [xopt, fopt, info]                  = self.subOptimizer.apply(x0);
%          [ok, msg, mv]                       = self.constraints.satisfied(xopt);
%          if ok
%             info.stopFlags.cstTolOk = true;
%             return
%          end
%          %
         info.stopFlags.maxFunEvalReached = false;
         xopt  = x0;
         fun   = self.initFun;
         %
         while self.nPen < self.maxPen && ~info.stopFlags.maxFunEvalReached && self.sigma < self.sigmaMax
            penSim                           = genetic.simulator(@(x) fun(x) + penF(x, self.sigma), self.subOptimizer.xDim, self.subOptimizer.nObj);
            penSim.nEval = self.subOptimizer.nEval;
%             penSim.maxFunEval = self.subOptimizer.maxFunEval;
%             opt                              = self.initOpt;
%             opt.maxFunEval                   = self.subOptimizer.maxFunEval - self.subOptimizer.nEval;
            self.subOptimizer                = genetic.optimizer.get(self.subOptimizer.methodName, penSim , self.initCst, self.initOpt);

            self.subOptimizer.constraints.skip  = self.toPenalise;

            [xopt, fopt, info]               = self.subOptimizer.apply(xopt);
            [ok,msg]                         = self.constraints.satisfied(xopt);
            fopt                             = fopt - penF(xopt, self.sigma);
            if ok
               info.stopFlags.cstTolOk = true;
               return
            else
               info.stopFlags.cstTolOk = false;
               info.vStopFlag          = [info.vStopFlag,' but ',msg];
               info.success            = false;
            end
            self.sigma                       = 2 * self.sigma;
         end
         %
         
      end
      %%
      function postEval(self, group)
      end
      %%
      function printReport(self, info)
         self.subOptimizer.printReport(info);
      end
   end
end
