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
classdef individual < handle
   % ===================================================================
   % ATTRIBUTES
   % ===================================================================
   properties (SetAccess=protected)
      % Current value and associated objective function
      value                   = [];
      previousValue           = [];
      objective               = [];
      previousObjective       = [];
      previousGrad            = [];
      grad                    = [];
      Hessian                 = [];
      % Memory
      best                    = [];
      previousBest            = [];
      % Controls
      nValueChange            = 0;
      nMoveWithoutImprove     = 0;
      bestObjectiveChange     = [];
      bestValueChange         = [];
      isAlive                 = true;
      neighbours              = [];
      %
      gradEstimator           = [];
      nImprove                = 0;
   end
   properties
      memorySize              = 1;
      id                      = 0;
      mother                  = [];
      data                    = [];
      nMaxValueChange         = inf;
      %
   end
   % ===================================================================
   % PUBLIC METHODS
   % ===================================================================
   methods
      %%
      % -------------------------------------------------------------------
      % INDIVIDUAL: GENERAL METHODS
      % -------------------------------------------------------------------      
      %% Constructor
      function obj = individual()
      end
      %%
      function xDim = xDim(self)
         xDim = length(self.value);
      end
      %% moveTo
      % moves the individual to a new value. Note that an individual can
      % move neither if it is dead nor if it is waiting for the previous
      % objective value. This has been made to avoid any coding error.
      function obj = moveTo(obj, newValue)
         if ~obj.isAlive
            error('Cannot move a dead individual')
         end
         %
         if ~isempty(obj.previousValue)
            obj.previousValue(:,2:end) = obj.previousValue(:,1:end-1);
            obj.previousValue(:,1)     = obj.value;
         else
            if ~isempty(obj.value)
               obj.previousValue       = zeros(size(obj.value,1), obj.memorySize);
               obj.previousValue(:,1)  = obj.value;
            end
         end
         obj.value               = newValue;
         obj.nValueChange        = obj.nValueChange + 1;
         if obj.nValueChange >= obj.nMaxValueChange
            obj.isAlive = false;
         end
      end
      %% getBestNeighbour
      % looks for the neighbour with the best objective value
      function x = getBestNeighbour(obj)
         x     = [];
         yBest = inf;
         for i = 1:length(obj.neighbours)
            neighbour = obj.mother.withID(obj.neighbours(i));
            if ~isempty(neighbour)
               if neighbour.objective < yBest
                  yBest = neighbour.objective;
                  x     = neighbour;
               end
            end
         end
      end
      %% tellObjective
      % tells the individual the value of the objective associated with the
      % current value of the individual and update its memory
      function obj = tellObjective(obj, newObjective, varargin)
%          obj.isWaitingObjective  = false;
         if ~isempty(obj.previousObjective)
               obj.previousObjective      = [obj.objective, obj.previousObjective(:,1:end -1)];
         else
            if ~isempty(obj.objective)
               obj.previousObjective      = zeros(size(obj.value,1), obj.memorySize);
               obj.previousObjective(:,1) = obj.objective;
            end
         end
         obj.objective           = newObjective;
         if length(varargin) >= 1
            obj.setGrad(varargin{1});
%             obj.grad    = varargin{1};
         end
         if length(varargin) >= 2
            obj.Hessian = varargin{2};
         end
         
         if ~isempty(obj.best)
            if obj.objective < obj.best.objective
               obj.nMoveWithoutImprove = 0;
               obj.nImprove      = obj.nImprove + 1;
               obj.previousBest  = obj.best;
               obj.store();
               if ~isempty(obj.previousBest)
                  obj.bestObjectiveChange = norm(obj.best.objective - obj.previousBest.objective,inf);
                  obj.bestValueChange     = norm(obj.best.value - obj.previousBest.value,inf);
               end
               % notify mother group that the individual has improved its
               % best value
               obj.mother.updateMemory(obj.best);
            else
               obj.nMoveWithoutImprove =  obj.nMoveWithoutImprove + 1;
            end
         else
            obj.nImprove   = obj.nImprove + 1;
            obj.store();
            % notify mother group that the individual has improved its
            % best value
            obj.mother.updateMemory(obj.best);
         end
         % indicates to the mother group that this individual is done
         % moving (it has its objective). It is used to increment the
         % number of generation.
%          obj.mother.individualHasMoved(obj.id);
         %
         if ~obj.isAlive
            obj.kill();
         end
      end
      %% store
      % stores the current value and objective of the individual in a
      % structure with the fields 'value', 'objective'.
      function store(obj)
         obj.best = struct('value',obj.value,'objective',obj.objective);
      end
      %% kill
      % kills the individual by removing it from its containing group.
      function kill(obj)
         obj.mother.killIndividuals(obj);
      end
      %% addNeighbour
      % adds a neighbour to the individual by its id.
      function addNeighbour(obj, id)
         obj.neighbours(end+1) = id;
      end
      %% setGrad
      function setGrad(self, g)
         self.previousGrad = self.grad;
         self.grad         = g;
      end
      %%
%       function reEvaluateBest(self, fun)
%          if ~isempty(self.best)
%             self.best.objective = fun(self.best.value);
%          end
%          if ~isempty(self.previousBest)
%             self.previousBest.objective = fun(self.previousBest.value);
%          end
%       end
   end
end