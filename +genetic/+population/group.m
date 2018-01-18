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
classdef group < handle
   % ===================================================================
   % ATTRIBUTES
   % ===================================================================
   properties (SetAccess = protected, Hidden = true)
      % Storage of individuals
      idMap                   = [];
      idxMap                  = [];
      nextIndex               = [];
      lastIndividualId        = 0;
      %
      trueChildIndexes        = [];
      nextChildIndex          = [];
   end
   properties (SetAccess = protected)
      mother                  = [];
      individuals
      len                     = 0;
      % child groups
      nChild                  = 0;
      child                   = [];
      % memory
      best                    = [];
      previousBest            = [];
      snapshotBest            = [];
      hasImproved             = false;
      % indicators for the group
      bestObjectiveChange     = [];
      bestValueChange         = [];
      nGen                    = 0;
      nGenLastChange          = 1;
      nImprove                = 0;
      freqImprove             = 0;
      %
      topologyStructure       = [];
      %
      individualsAreEncoded   = false;
   end
   properties
      discardUpdate = false;
      overlay                 = [];
      name                    = '';
      preventUpdate           = false;
      %
      id                      = [];
      % user data
      data                    = struct();   
   end
   %%
   % ======================================================================
   % PUBLIC METHODS
   % ======================================================================
   methods
      %%
      % -------------------------------------------------------------------
      % GROUP: GENERAL METHODS
      % -------------------------------------------------------------------
      %% Constructor
      function obj = group(initGroupSize)
         persistent id
         if ~exist('id','var')
            id = 1;
         else
            id = id+1;
         end
         if nargin == 1
            obj.spawnIndividual(initGroupSize);
         end
         obj.id = id;
      end
      %%
      function xDim = xDim(self)
         xDim = [];
         if self.len > 0
            xDim = self.at(1).xDim;
         end
      end
      %%
      function yDim = yDim(self)
         yDim = size(self.at(1).objective,1);
      end
      %% hasImproved
      function hasImproved_(self, replace)
         if nargin < 2
            replace = false;
         end
         if ~isempty(self.snapshotBest)
            out = self.best.objective < self.snapshotBest.objective ;
         else
            out = true;
         end
         if replace
            self.snapshotBest = self.best;
         end
         
         self.hasImproved = out;
      end
      %% initTopology
      % initialises the matrix representing the (initial) topology of the
      % group.
      function initTopology(obj, topologyType, assignTopo)
         if obj.len < 2
            error('There must be at least 2 particles to initialise the topology of the group');
         end
         if nargin < 2
            topologyType = obj.topologyType;
         end
         if nargin < 3
            assignTopo = false;
         end
         %
         switch topologyType
            case 'full'
               T = genetic.population.topology.full(obj.len);
            case 'vonn'
               T = genetic.population.topology.vonNeumann(obj.len);
            case 'ring'
               T = genetic.population.topology.ring(obj.len);
         end
         obj.topologyStructure = T;

         if assignTopo
            obj.assignTopology()
         end
      end
      %% assignTopology
      % assign the neighbours to each individual based on the topology
      % matrix.
      function assignTopology(obj)
         for i = 1:obj.len
            xi                = obj.at(i);
            neighboursNumber  = obj.topologyStructure{i};
            for j = 1:length(neighboursNumber)
               xjId = obj.at(neighboursNumber(j)).id;
               xi.addNeighbour(xjId);
            end
         end
      end
      %% updateMemory
      function updateMemory(obj, candidate)
         % compares a given structure with the current best structure and
         % updates the latter if needed. If a mother group exists, update its
         % memory too.
         if ~isempty(obj.best)
            if candidate.objective < obj.best.objective
               obj.nImprove         = obj.nImprove + 1;
               obj.previousBest     = obj.best;
               obj.best             = candidate;
               if ~isempty(obj.mother) && ~obj.preventUpdate
                  obj.mother.updateMemory(obj.best);
               end
               if ~isempty(obj.previousBest)
                  obj.bestObjectiveChange    = norm(obj.best.objective - obj.previousBest.objective,inf);
                  obj.bestValueChange        = norm(obj.best.value - obj.previousBest.value,inf);
               end
            end
         else
            obj.nImprove      = obj.nImprove + 1;
            obj.best          = candidate;
            if ~isempty(obj.mother) && ~obj.preventUpdate
               obj.mother.updateMemory(obj.best);
            end
         end
      end
      %% applyToIndividuals
      function out = applyToIndividuals(obj, call, indexes, asMat)
         % applies some function to all the individuals of the group and
         % returns the result in a matrix (default) or a cell.
         % This routine is for example used when getting the values or 
         % objective functions of the population.
         if nargin < 3 || isempty(indexes)
            indexes  = 1:obj.len;
            asMat    = true;
         end
         if nargin < 4 || isempty(asMat)
            asMat = true;
         end
         if ~iscell(call)
            call = {call};
         end
         nCalls   = length(call);
         out      = cell(1,nCalls);
         for i = 1:length(indexes)
            ind = obj.at(indexes(i));
            for j = 1:nCalls
               out{j}{i} = call{j}(ind);
            end
         end
         if nCalls == 1
            out = out{1};
            if asMat
               out = cell2mat(out);
            end
         else
            if asMat
               for i = 1:nCalls
                  out{i} = cell2mat(out{i});
               end
            end
         end
      end
      %% getMinMax
      function [m,M] = getMinMax(self)
         % returns the minimal and maximal values of the objective
         % functions of the current population. 
         % For multi-objective functions, the min-max values are considered
         % over each objective.
         Y     = self.getObjective();
         [m,M] = genetic.population.group.getExtremum(Y);
      end
      %% getObjective
      function out = getObjective(obj, varargin)
         % returns the values of the objective functions of the current
         % population.
         extract  = @(I) I.objective;
         out      = obj.applyToIndividuals(extract, varargin{:});
      end
      %% getValue
      function out = getValue(obj, varargin)
         % returns the value of all the individuals of the group, as a cell
         % (default) or as a matrix
         extract  = @(I) I.value;
         out      = obj.applyToIndividuals(extract,varargin{:});
      end
      %% getBestValues
      function out = getBestValues(obj, varargin)
         extract  = @(I) I.best.value;
         out      = obj.applyToIndividuals(extract,varargin{:});
      end
      %% getBestNValues
      function out = getBestNValues(obj, varargin)
         % returns the values of the best neighbours of each individual
         extract  = @(I) I.getBestNeighbour().value;
         out      = obj.applyToIndividuals(extract,varargin{:});
      end
      %% getBestObjectives
      function out = getBestObjectives(obj, varargin)
         extract  = @(I) I.best.objective;
         out      = obj.applyToIndividuals(extract,varargin{:});
      end
      %% getBestNObjectives
      function out = getBestNObjectives(obj, varargin)
         % returns the objectives of the best neighbours of each individual
         extract  = @(I) I.getBestNeighbour().objective;
         out      = obj.applyToIndividuals(extract,varargin{:});
      end
      %% getGrads
      function out = getGrads(obj, varargin)
         extract  = @(I) I.grad;
         out      = obj.applyToIndividuals(extract, varargin{:});
      end
      %%
      function [xopt, fopt] = getAbsoluteBest(self)
         xopt = [];
         fopt = [];
         if ~isempty(self.mother)
            [xopt, fopt] = self.mother.getAbsoluteBest();
         else
            if ~isempty(self.best)
               xopt = self.best.value;
               fopt = self.best.objective;
%             else
%                xopt = self.getValue();
%                fopt = self.getObjective();
            end
         end
      end
      %%
      function g = getAncestor(self)
         if ~isempty(self.mother)
            g = self.mother.getAncestor;
         else
            g = self;
         end
      end
      %% sort
      function S = sort(self, objective)
         if nargin < 2 
            Y = self.getObjective;
            N = self.len;
         else
            Y = objective;
            N = size(Y,2);
         end
         [Ys, Is] = sort(Y,2,'ascend');
         S        = cell(N,1);
         lastBest = inf;
         for i = 1:N
            if Ys(i) == lastBest
               S{end}   = [S{end},Is(i)];
            else
               S{end+1} = Is(i);
            end
            lastBest = Ys(i);
         end
         S = genetic.tools.removeEmptyCell(S);
      end
      %% moveTo
      function moveTo(obj, val)
         % moves the individuals of the group to a new position
         for i = 1:obj.len
            obj.at(i).moveTo(val(:,i));
         end
      end
      %% tellObjective
      function tellObjective(self,y, varargin)
         for i = 1:self.len
            self.at(i).tellObjective(y(:,i), varargin{:});
         end
         self.hasImproved_(true); % Test if improvement + save snapshot
      end
      %% kill
      % kills the group and removes it from the child list of its mother
      % (if any).
      function kill(obj)
         if ~isempty(obj.mother)
            obj.mother.killChildGroup(obj);
         end
      end
      %%
      % -------------------------------------------------------------------
      % HANDLING OF INDIVIDUALS
      % -------------------------------------------------------------------
      %% spawnIndividual
      % spawns a given number (default 1) of individuals and adds them to
      % the group
      function spawnIndividual(obj, nIndividuals, values)
         if nargin < 2
            nIndividuals = 1;
            values = [];
         end
         if nargin < 3
            values = [];
         end
         for i = 1:nIndividuals
            x = genetic.population.individual();
            if ~isempty(values)
               x.moveTo(values(:,i));
            end
            obj.addIndividual(x);
         end
      end
      %% addIndividual
      % adds a given individual to the group.
      function addIndividual(obj, ind)
         ID                               = obj.getNextIndividualId();
         idx                              = obj.getNextIndividualPosition();
         obj.addMapping(ID,idx);
         ind.id                           = ID;
         ind.mother                       = obj;
         obj.individuals{idx}             = ind;
         if ~isempty(ind.objective)
            obj.updateMemory(ind.best);
         end
         obj.len                          = obj.len + 1;
      end
      %% at
      % extracts the i-th individual from the group.
      % As the individuals can be removed and added, the id of the
      % extracted individual does not necessarily corresponds to its
      % position <i>.
      function individual = at(obj, i)
         if i <= obj.len
            idx         = obj.idxMap(i);
            individual  = obj.individuals{idx};
         else
            error('Individual does not exist')
         end
      end
      %% withID
      % extracts an individual given its id
      function individual = withID(obj, id)
         idx = obj.idxMap(obj.idMap == id);
         if ~isempty(idx)
            individual = obj.individuals{idx};
         else
            error('Individual not found');
         end
      end
      %% killIndividual
      % removes a given individual from the group.
      function killIndividuals(obj, X)
         if isa(X, 'genetic.population.individual')
            X = {X};
         end
         if isa(X, 'double')
            % If the input is a set of indexes, then each corresponding 
            % individual is killed.
            to_kill = sort(X,'descend');
            for i = to_kill
               obj.at(i).kill();
            end
         elseif isa(X, 'cell')
            % If the input is a set of individuals, then each individual is
            % killed.
            for i = 1:length(X)
               % /!\ Ajouter un check pour être sûr que x appartient bien à ce
               % groupe
               x                    = X{i};
               obj.len              = obj.len -1;
               [idx, idxMapping]    = obj.getIndividualPosition(x.id);
               obj.individuals{idx} = [];
               obj.nextIndex(end+1) = idx;
               obj.deleteMapping(idxMapping);
%                obj.nThatMustMove    = obj.nThatMustMove - 1;
            end
         end
      end
      %%
      % -------------------------------------------------------------------
      % HANDLING OF CHILD GROUPS
      % -------------------------------------------------------------------
      %% spawnGroup
      % spawns a new group with a given optimizer and given individuals
      % (both are optional).
      function newGroup = spawnGroup(obj, optimizer, ind)
         idx                  = obj.getNextChildIdx();
         obj.addChildIndex(idx);
         newGroup             = genetic.population.group();
         newGroup.id          = idx;
         newGroup.individualsAreEncoded  = obj.individualsAreEncoded;
%          newGroup.evalFun     = obj.evalFun;
%          newGroup.xDim        = obj.xDim;
%          newGroup.minValue    = obj.minValue;
%          newGroup.maxValue    = obj.maxValue;
% %          newGroup.maxFunEval  = obj.maxFunEval;
         newGroup.data        = obj.data;
         newGroup.mother      = obj;
         if nargin >= 2
            newGroup.optimizer   = optimizer;
            if nargin == 3
               if isa(ind,'cell')
                  for i = 1:length(ind)
                     x = ind{i}.copy();
                     newGroup.addIndividual(x);
                  end
               else
                  x = ind.copy();
                  newGroup.addIndividual(x);
               end
            end
         end
         obj.child{idx} = newGroup;
         obj.nChild     = obj.nChild + 1;
      end
      %% addChild
      % adds an existing group as a child of another group.
      function addChild(obj, group)
         idx                  = obj.getNextChildIdx();
         obj.addChildIndex(idx);
         obj.child{idx}       = group;
         group.setMother(obj);
         group.id             = idx;
         obj.nChild           = obj.nChild + 1;
      end
      %% setMother
      % sets the mother group.
      function setMother(obj, mother)
         obj.mother = mother;
      end
      %% reEvaluateBest
      % re-evaluates the solution stored in the memory and erase the stored
      % value of the objective with the result. This is used when the
      % objective function is update through the optimisation process (for
      % instance with penalisation) and enables to make the past solution
      % comparable to the new ones.
      function reEvaluateBest(self, fun)
         if ~isempty(self.best)
            self.best.objective = fun(self.best.value);
         end
         if ~isempty(self.previousBest)
            self.previousBest.objective = fun(self.previousBest.value);
         end
      end
      %% gAt
      % returns the i-th created child group
      function childGroups = gAt(obj, idx)
         if max(idx)>obj.nChild
            error('child group does not exist');
         end
         ID          = obj.trueChildIndexes(idx);
         n           = length(idx);
         childGroups = cell(n,1);
         for i = 1:n
            childGroups{i} = obj.child{ID(i)};
         end
         if n == 1
            childGroups = childGroups{1};
         end
      end
      %% killChildGroup
      % remove a given child group from the group.
      % /!\ vérifier que le groupe fourni est bien un enfant du groupe
      % courant.
      function killChildGroup(obj, g)
         obj.child{g.id}            = [];
         obj.deleteChildIndex(g.id);
         obj.nextChildIndex(end+1)  = g.id;
         obj.nChild                 = obj.nChild - 1;
      end
   end % end of public methods
   %%
   % ======================================================================
   % PRIVATE METHODS
   % ======================================================================
   methods (Access=private)
      %%
      % -------------------------------------------------------------------
      % GENERAL TOOLS
      % -------------------------------------------------------------------
      %% structureToTopology
      % converts a generic topology into a topology based on the IDs of the
      % particles in the group.
      function newT = structureToTopology(obj, T)
         newT  = T;
         for i = 1:obj.len
            xiNeighbours = T{i};
            for j = 1:length(xiNeighbours)
               newT{i}(j) = obj.at(xiNeighbours(j)).id;
            end
         end
      end
      %%
      % -------------------------------------------------------------------
      % INDIVIDUALS
      % -------------------------------------------------------------------
      %% getIndividualPosition
      % get the position of an individual in the individuals list from its
      % ID.
      function [idx,idxMapping] = getIndividualPosition(obj, id)
         idxMapping  = find(obj.idMap == id);
         idx         = obj.idxMap(idxMapping);
      end
      %% addMapping
      % store the mapping between the ID of an individual and its position
      % idx in the individuals list.
      function addMapping(obj, ID, idx)
         obj.idMap(end + 1)   = ID;
         obj.idxMap(end + 1)  = idx;
      end
      %% deleteMapping
      % removes an index from the list of non-empty individuals. This
      % routine is called when an individual is deleted.
      function deleteMapping(obj,idx)
         obj.idxMap(idx) = [];
         obj.idMap(idx)  = [];
      end
      %% getNextIndividualId
      % generates a new id for a new individual. Note that it does not
      % correspond to the position in the individuals list.
      function id = getNextIndividualId(obj)
         id                   = obj.lastIndividualId + 1;
         obj.lastIndividualId = id;
      end
      %% getNextIndividualPosition
      % determines the position at which a new individual must be stored in
      % the individuals list.
      function idx = getNextIndividualPosition(obj)
         if ~isempty(obj.nextIndex)
            idx               = obj.nextIndex(1);
            obj.nextIndex(1)  = [];
         else
            idx               = obj.len + 1;
         end
      end
      %%
      % -------------------------------------------------------------------
      % CHILD GROUPS
      % -------------------------------------------------------------------
      %% getNextChildIdx
      % determines the position at which a new child group must be stored.
      function idx = getNextChildIdx(obj)
         if ~isempty(obj.nextChildIndex)
            idx                     = obj.nextChildIndex(1);
            obj.nextChildIndex(1)   = [];
         else
            idx                     = obj.nChild + 1;
         end
      end
      %% addGroupIndex
      % adds an index to the list of non-empty child groups. This routine
      % is called when a child group is created
      function addChildIndex(obj,id)
         obj.trueChildIndexes(end+1) = id;
      end
      %% deleteChildIndex
      % removes an index from the list of non-empty child groups. This
      % routine is called when a child group is deleted.
      function deleteChildIndex(obj,id)
         obj.trueChildIndexes(obj.trueChildIndexes == id) = [];
      end
   end % end of private methods
   % ======================================================================
   % STATIC METHODS
   % ======================================================================
   methods (Static)
      %%
      function [m,M] = getExtremum(Y)
         m = min(Y,[],2);
         M = max(Y,[],2);
      end
      %%
      function x = initValue(initMethod, popSize, x0, xMin, xMax, regularize, defaultInf)
         if nargin < 7 
            defaultInf = 1e3;
         end
         % Replace infinite value
         xMin(xMin == -inf)   = -defaultInf;
         xMax(xMax == inf)    = defaultInf;
         %
         xDim = length(xMin);
         if  popSize > 0
            switch initMethod
               case 'rnd'
                  x = genetic.population.init.random( popSize, xMin, xMax, xDim);
               case 'rem'
                  x = genetic.population.init.randomRem( popSize, xMax, xMax, xDim);
               case 'ham'
                  x = genetic.population.init.hammersley( popSize, xMin, xMax, xDim);
               case 'hal'
                  x = genetic.population.init.halton( popSize, xMax, xMax, xDim);
               case 'vor'
                  x = genetic.population.init.voronoi( popSize, xMin, xMax, xDim);
            end
            if regularize
               x = genetic.population.init.regularise(x,  popSize, xMin, xMax, xDim);
            end
         end
         if ~isempty(x0)
            m        = min(size(x0,2),popSize);            
            x(:,1:m) = x0;
         end
      end
      %% select
      function [X,idx] = select(n, X, Y, method, varargin)
         if nargin < 4 
            method = 'rand';
         end
         nX = size(X,2);
         switch method
            case 'rand'
               idx   = randperm(nX,n);               
            otherwise
               error('Not implemented')
         end
         X = X(:,idx);
      end
      %% 
      function X = tabulate(i)
         counts   = histcounts(i,'BinMethod','integers');
         total    = sum(counts);
         X        = [(1:max(i))', counts', 100*counts'./total];
         
      end
     % getFreeIdx
%       function idx = getFreeIdx(c)
%          % finds an empty space in the list of individuals and returns the
%          % corresponding index. If no space is available
%          if isempty(c)
%             idx = 1;
%          else
%             idx = find(cellfun('isempty',c),1);
%             if isempty(idx)
%                idx = length(c) + 1;
%             end
%          end
%       end
   end
end
