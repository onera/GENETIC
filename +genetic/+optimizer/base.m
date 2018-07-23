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
classdef base < handle
   %
   properties
      simulator         = [];
      %
      methodName        = '';
      longName          = '';
      verbosity         = 0;
      %
      maxFunEval        = 5000;
      maxGen            = inf; 
      % Constraints
      constraints       = [];
      defaultInfinity   = 1e3;
      needFiniteBounds  = false;
      handledConstraints = {};
      %
      stopTests         = {struct('name','maxFunEvalReached', 'cleanStop', false, 'condByConstraints', false), struct('name','maxGenReached','cleanStop',false,'condByConstraints',false)};
      stopFlags         = [];
      vStopFlag         = '';
      success           = false;
      % Population
      popSize           = [];
      initMethod        = 'rnd';
      regulariseInit    = false;
      %
      toDisplay         = {};
      iterFlag          = '';
      %
      printLvl          = 1;
%       modPrintIter      = 1;
      modPrintHead      = 50;
      printOnlyImp      = false;
      printCtr = 0;
   end

   methods
      %% Constructor
      function self = base(simulator, constraints, options)
         if nargin < 1
            error('At least the simulator must be provided')
         end
         self.simulator    = simulator;
         self.constraints  = constraints;
         %
         if nargin == 3
            self.fillAttributes(options);
         end
         %
         self.constraints.adjustBoundsDimension(self.xDim);
         %
         self.simulator.maxFunEval = self.maxFunEval;
         if isempty(self.popSize)
            self.setDefaultPopSize();
         end
         %
         if genetic.tools.params('verbosity') == 2
            self.printOnlyImp = true;
         end
      end
      %%
      function optimizer = wrapper4constraints(self, constraints, options)
         cstNotHandled = setdiff(self.constraints.list(), self.handledConstraints);
         if ~isempty(cstNotHandled)
            optimizer = genetic.optimizer.penalty(self, cstNotHandled, constraints, options);
         else
            optimizer = self;
         end
      end
      %% setOptions
      function fillAttributes(self, options, quiet)
         if nargin < 3
            quiet = false;
         end
         if ~isempty(options)
            fieldsNames = fieldnames(options);
            for i = 1:length(fieldsNames)
               try
                  self.(fieldsNames{i}) = options.(fieldsNames{i});
               catch
                  if ~quiet
                     fprintf('''%s'' is not an attribute of %s, discarding.\n',fieldsNames{i},class(self));
                  end
               end
            end
         end
      end
      % ===================================================================
      %                       SIMULATOR-RELATED
      % ===================================================================
      %% nEval getter
      function nEval = nEval(self)
         nEval = self.simulator.nEval;
      end
      %% xDim getter
      function xDim = xDim(self)
         xDim = self.simulator.xDim;
      end
      %% nObj getter
      function nObj = nObj(self)
         nObj = self.simulator.nobj;
      end
      %% xMin getter
      function xMin = xMin(self)
         xMin = self.constraints.xMin;
         if self.needFiniteBounds
            xMin(xMin == -inf) = -self.defaultInfinity;
         end
      end
      %% xMax getter
      function xMax = xMax(self)
         xMax = self.constraints.xMax;
         if self.needFiniteBounds
            xMax(xMax == inf) = self.defaultInfinity;
         end
      end
      %% isBoxBounded
      function out = isBoxBounded(self)
         out = self.constraints.containsBox();
      end
      %% isUnconstrained
      function out = isUnconstrained(self)
         out = self.constraints.empty();
      end
      % ===================================================================
      %                          DISPLAY STUFF
      % ===================================================================
      %% addDisplayElement
      function addDisplayElement(self,varargin)
         for i = 1:length(varargin)
            self.toDisplay{end+1} = varargin{i};
         end
      end
      %%
      function removeDisplay(self, key)
         for i = 1:length(self.toDisplay)
            if strcmp(self.toDisplay{i}.name,key)
               self.toDisplay{i} = [];
            end
         end
         self.toDisplay = genetic.tools.removeEmptyCell(self.toDisplay);
      end
      %% printImp
      function out = printImp(self, group)
         out = '';
         if group.hasImproved
            out = '*';
         end
      end
      %% printnGen
      function out = printnGen(self, group)
         out = sprintf('%d',self.nGen);
      end
      %% printNeval
      function out = printNeval(self, group)
         out = sprintf('%d/%d',  self.nEval, self.maxFunEval);
      end      
      %% printHeader
      function printHeader(self, x0Provided)
         if genetic.tools.params('verbosity')< self.printLvl
            return
         end
         al = ['r','c','l'];
         d  = [16;1;100];
         %
%          genetic.tools.print(1, ' Genetic Toolbox');
         %
         el = {'method', sprintf('%s (key: ''%s'')',self.longName, self.methodName)};
         genetic.tools.print(1,genetic.tools.printCol({el{1},':',el{2}}, al, d));
         %
         el = {'dimensions', sprintf('x:%d, y:%d',self.xDim, self.nObj)};
         genetic.tools.print(1,genetic.tools.printCol({el{1},':',el{2}}, al, d));
         %
         if x0Provided
            el = {'init', 'user-defined'};
         else
            el = {'init', 'random'};
         end
         genetic.tools.print(1,genetic.tools.printCol({el{1},':',el{2}}, al, d));
      end
      %% printReport
      function printReport(self, info)
         if genetic.tools.params('verbosity')< self.printLvl
            return
         end
         al = ['r','c','l'];
         d  = [16;1;100];
         %
         el = {'success', sprintf('%d',info.success)};
         genetic.tools.print(1,genetic.tools.printCol({el{1},':',el{2}}, al, d));
         %
         el = {'nEval', sprintf('%d',info.nEval)};
         genetic.tools.print(1,genetic.tools.printCol({el{1},':',el{2}}, al, d));
         %
         el = {'flag', info.vStopFlag};
         genetic.tools.print(1,genetic.tools.printCol({el{1},':',el{2}}, al, d));
         %
         el = {'time', sprintf('%.1fs', info.elapsedTime)};
         genetic.tools.print(1,genetic.tools.printCol({el{1},':',el{2}}, al, d));

      end
      %% printStart
      function printStart(self)
         if genetic.tools.params('verbosity')< self.printLvl
            return
         end
         genetic.tools.print(1, ' Optimisation started...');
      end
      %% printEnd
      function printEnd(self)
         if genetic.tools.params('verbosity')< self.printLvl
            return
         end
         if self.success
            genetic.tools.print(1, ' Optimisation successful.');
         else
            genetic.tools.print(1, ' Optimisation done.');
         end
      end
      %% printIterHead
      function printIterHead(self)
         if genetic.tools.params('verbosity') < self.printLvl + 1
            return
         end
         ntd   = length(self.toDisplay);
         el    = cell(ntd,1);
         al    = '';
         d     = zeros(ntd,1);
         for i = 1:ntd
            el{i} = self.toDisplay{i}.name;
            al    = [al;self.toDisplay{i}.align];
            d(i)  = self.toDisplay{i}.dim;
         end
         if ntd > 0
            headerString   = genetic.tools.printCol(el, al, d);
            genetic.tools.print(2, ['\n',headerString]);
         end
      end
      %% printIter
      function printIter(self, group)
         if genetic.tools.params('verbosity') < self.printLvl + 1
            return
         end
         if self.printOnlyImp && ~group.hasImproved
            return
         end
         self.printCtr  = self.printCtr + 1;
         ntd            = length(self.toDisplay);
         el             = cell(ntd,1);
         al             = '';
         d              = zeros(ntd,1);
         for i = 1:ntd
            el{i} = self.(self.toDisplay{i}.call)(group);
            al    = [al;self.toDisplay{i}.align];
            d(i)  = self.toDisplay{i}.dim;
         end
         if ntd > 0
            iterString = genetic.tools.printCol(el, al, d);
            genetic.tools.print(2, iterString);
         end
      end
      % ===================================================================
      %                          STOPPING CRITERIA
      % ===================================================================
      %% addStopTest
      function addStopTest(self, varargin)
         for i = 1:length(varargin)
            self.stopTests{end+1} = varargin{i};
         end
      end
      %% stop
      % Tests a list of stopping criteria
      function out = stop(self, group)
         out            = false;
         checkCst       = true;
         %
         self.stopFlags = struct();
         self.vStopFlag = '';
         for i = 1:length(self.stopTests)
            t                       = self.stopTests{i};
            [stop, msg]             = self.(t.name)(group); % overhead lent
            self.stopFlags.(t.name) = stop;
            self.vStopFlag          = [self.vStopFlag, msg];
            if stop && t.condByConstraints
               if checkCst
                  checkCst       = false;
                  [cstOk,msgC]   = self.cstSatisfied(group);
                  self.vStopFlag = [self.vStopFlag, msgC];
               end
               stop = stop & cstOk;
            end
            out = out | stop;
            if stop && t.cleanStop
               self.success = true;
            end
         end
      end
      %% cstSatisfied
      function [out, msg] = cstSatisfied(self, group)
         out = false;
         msg = '';
         if ~isempty(group.best.value)
            x           = group.best.value;
            [out,msg]   = self.constraints.satisfied(x);
         end
      end
      %% maxEvalReached
      % tests whether the maximum number of evaluations has been reached
      function [out, msg] = maxFunEvalReached(self, dummy)
         [out, msg] = self.simulator.maxFunEvalReached();
      end
      %% maxGenReached
      function [out, msg] = maxGenReached(self, dummy)
         out = self.nGen >= self.maxGen;
         msg = '';
         if out
            msg = sprintf('Maximum number of generation reached: %d',self.nGen);
         end
      end
      % ===================================================================
      %                          MISC
      % ===================================================================
      function x0 = initPopulation(self, x0)
         x0 = genetic.population.group.initValue(self.initMethod, self.popSize, x0, self.xMin, self.xMax, self.regulariseInit, self.defaultInfinity);
      end
      %%
      function info = fillInfo(self, group)
         info              = struct();
         info.success      = self.success;
         info.nEval        = self.nEval;
         info.stopFlags    = self.stopFlags;
         info.vStopFlag    = self.vStopFlag;
         info.nImprove     = group.nImprove;
         info.freqImprove  = info.nImprove/info.nEval;
         info.lastPop      = struct('x', group.getValue(), 'y', group.getObjective());
      end
      %%
      function setDefaultPopSize(self)
         self.popSize = round(10 + 2 *sqrt(self.xDim));
      end
   end
end
