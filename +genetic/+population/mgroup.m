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
classdef mgroup < genetic.population.group
   properties
      maxFrontLen          = inf;
      orderRules           = 'pareto';
      idealPoint           = [];
      nadirPoint           = [];
      metric               = '';
      metricValue          = Inf;
      previousFront        = [];
   end
   methods
      %% Constructor
      function obj = mgroup(initGroupSize)
         if nargin == 1
            gs = initGroupSize;
         else
            gs = 0;
         end
         obj@genetic.population.group(gs);
      end
      %% hasImproved
      function hasImproved_(self, replace)
         if nargin < 2
            replace = false;
         end
         if ~isempty(self.snapshotBest)
            out = size(self.best.objective,2) ~= size(self.snapshotBest.objective,2) || any(any(self.best.objective ~= self.snapshotBest.objective));
         else
            out = true;
         end
         if replace
            self.snapshotBest = self.best;
         end
         
         self.hasImproved = out;
      end
      %% updateMemory
      % compares a given structure with the current best structure and
      % updates the latter if needed. If a mother group exists, update its
      % memory too.
      function updateMemory(obj, candidate)
%          if obj.discardUpdate
%             return
%          end
         if isempty(obj.mother)
            % if the group has no mother group, then its memory is the
            % front. The candidate is tested with the latter.
            if ~isempty(obj.best)
               % if the front is not empty, the candidate must be compared
               % with the front to see if it fits
               front             = obj.best;
               lenFront          = size(front.objective,2);
               %
               candidateObj      = repmat(candidate.objective,[1,lenFront]);
               % Test first if the point is not already in the front
               alreadyInFront    = any(all(candidateObj == front.objective,1),1);
               if alreadyInFront
                  return
               end 
               %
               switch obj.orderRules
                  case 'pareto'
                     candidateDominates                  = all(candidateObj <= front.objective,1) & any(candidateObj ~= front.objective,1);
                     candidateIsDominated                = all(front.objective <= candidateObj,1) & any(candidateObj ~= front.objective,1);
               end
               %
               %
               if any(candidateDominates) || ~any(candidateIsDominated)
                  front.objective(:,candidateDominates)  = [];
                  front.value(:,candidateDominates)      = [];
                  front.objective(:,end+1)               = candidate.objective;
                  front.value(:,end+1)                   = candidate.value;
                  obj.nImprove                           = obj.nImprove + 1;
                  obj.hasImproved                        = true;
                  obj.previousBest                       = obj.best;
                  obj.best                               = front;
               end
               %
            else
               % if the front is empty, the candidate becomes the front
               obj.best          = candidate;
               obj.nImprove      = obj.nImprove + 1;
               obj.hasImproved   = true;
            end
            if size(obj.best.objective,2) > obj.maxFrontLen
               obj.decimateFront();
            end
         else
            % if the group has a mother group, then the memory is in the
            % mother group (or above). The candidate is just passed.
            obj.mother.updateMemory(candidate);
         end
      end
      %% updateMetric
      function updateMetric(obj)
         % Isolate mother group
         group             = obj.topMother();
         % Do not calculate the metric value if there is no best field in
         % group
         if isempty(group.best)
            return;
         end
         %
         % Otherwise: 1/ stock the best objective values at current
         %               iteration and the previous ones (if they exist);
         %            2/ and update both ideal and nadir points of the
         %               group.
         yBest             = group.best.objective;
         %
         if ~isempty(group.previousFront)
            yPrevBest = group.previousFront.objective;
         else
            yPrevBest = [];
         end
         %
         group.idealPoint  = min([min(group.best.objective,[],2) group.idealPoint],[],2);
         group.nadirPoint  = max([max(group.best.objective,[],2) group.nadirPoint],[],2);
         %
         idealPt           = group.idealPoint;
         nadirPt           = group.nadirPoint;
         % Switch between the available indicators
         switch group.metric
            case 'HVM'
               val = genetic.population.metrics.hyperVolume(yBest, nadirPt);
            case 'DSM'
               val = genetic.population.metrics.distribMetric(yBest);
            case 'EXT'
               val = genetic.population.metrics.extentMetric(yBest);
            case 'SPA'
               val = genetic.population.metrics.spacingMetric(yBest);
            case 'MCD'
               val = genetic.tools.metrics.maxCrowdingDistance(yBest);
            case 'DVM'
               val = genetic.tools.metrics.diversityMetric(yBest);
            case 'NDC'
               val = genetic.tools.metrics.nDistinctChoices(yBest);
            case 'MCL'
               val = genetic.tools.metrics.muClusterMetric(yBest);
            case 'CRM'
               val = genetic.tools.metrics.coverRateMetric(yBest);
            case 'MSM'
               val = genetic.tools.metrics.maxSpreadMetric(yBest);
            case 'MSU'
               val = genetic.tools.metrics.minSumMetric(yBest);
            case 'SUM'
               val = genetic.tools.metrics.sumMinMetric(yBest);
            case 'SRM'
               val = genetic.tools.metrics.sumRangeMetric(yBest);
            case 'BSS'
               val = genetic.tools.metrics.bsSpreadMetric(yBest);
            case 'MDG'
               val = genetic.tools.metrics.minDistanceGraph(yBest);
            case 'FSM'
               val = genetic.tools.metrics.frontSpreadMetric(yBest);
            case 'SPM'
               val = genetic.tools.metrics.spreadMeasure(yBest);
            case 'ENT'
               if ~isempty(idealPt) && ~isempty(nadirPt)
                  val = genetic.tools.metrics.entropyMetric(yBest, idealPt, nadirPt);
               else
                  val = Inf;
               end
            case 'OPS'
               if ~isempty(idealPt) && ~isempty(nadirPt)
                  val = genetic.tools.metrics.overallParetoSpread(yBest, idealPt, nadirPt);
               else
                  val = Inf;
               end
            case 'RNI'
               Y     = group.getObjective();
               F     = genetic.population.tools.getNonDominatedIndividuals(Y);
               val   = genetic.tools.metrics.ratioNonDominatedIndividuals(Y, F);
            case 'DQM'
               xMin  = group.optimizer.xMin;
               xMax  = group.optimizer.xMax;
               val   = genetic.tools.metrics.dQuantifierMetric(yBest, xMin, xMax);
            case 'CSR'
               if ~isempty(yPrevBest)
                  val = genetic.tools.metrics.consolidationRatio(yBest, yPrevBest);
               else
                  val = Inf;
               end
            case 'CTR'
               if ~isempty(yPrevBest)
                  val = genetic.tools.metrics.contributionRatio(yBest, yPrevBest);
               else
                  val = Inf;
               end
            case 'EPS'
               if ~isempty(yPrevBest)
                  val = genetic.tools.metrics.epsilonIndicator(yBest, yPrevBest);
               else
                  val = Inf;
               end
            case 'SCM'
               if ~isempty(yPrevBest)
                  val = genetic.tools.metrics.setConvergenceMetric(yBest, yPrevBest);
               else
                  val = Inf;
               end
            case 'MDR'
               if ~isempty(yPrevBest)
                  val = genetic.tools.metrics.mutualDominationRate(yBest, yPrevBest);
               else
                  val = Inf;
               end
            case 'TSC'
               if ~isempty(yPrevBest)
                  val = genetic.tools.metrics.twoSetCoverage(yBest, yPrevBest);
               else
                  val = Inf;
               end
            case 'VMI'
               if ~isempty(yPrevBest)
                  val = genetic.tools.metrics.volumeMeasureIndicator(yBest, yPrevBest);
               else
                  val = Inf;
               end
            otherwise
               val = Inf;
         end
         %
         group.metricValue   = val;
         group.previousFront = group.best;
         %
      end
      %% sort
      % /!\ This probably involves unnecessary comparisons
      function [F,Ys] = sort(self)
         Y        = self.getObjective();
         [F, Ys]  = genetic.population.mgroup.paretoSort(Y);
      end      
      %% decimateFront
      % /!\ experimental
      function decimateFront(obj)
         if size(obj.best.objective, 2) - obj.maxFrontLen > 0
            d                    = genetic.optimizer.mu.nsga2.crowdingDist(obj.best.objective);
            [~, idx]             = sort(d,'descend');
            obj.best.objective   = obj.best.objective(:,idx(1:obj.maxFrontLen));
            obj.best.value       = obj.best.value(:,idx(1:obj.maxFrontLen));
         end  
      end
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
            x = genetic.population.mindividual();
            if ~isempty(values)
               x.moveTo(values(:,i));
            end
            obj.addIndividual(x);
         end         
      end
      %% topMother
      function m = topMother(self)
         if isempty(self.mother)
            m = self;
         else
            m = self.mother.topMother;
         end
      end
   end
   %%
   % ===================================================================
   % STATIC METHODS
   % ===================================================================   
   methods (Static)
      %% getNonDominatedIndividuals
      function F = getNonDominatedIndividuals(Y)
         n              = size(Y,2);
         candidatesMask = ones(n,1);
         % candidatesIdx  = 1:n;
         F              = [];
         iter           = 0;
         while any(candidatesMask~=0)
            iter                                = iter + 1;
            i                                   = find(candidatesMask,1);
            candidate                           = Y(:,i);
            %
            C                                   = repmat(candidate,1,n);
            diffTest                            = any(C ~= Y,1);
            candidateDominates                  = all(C <= Y,1) & diffTest;
            %
            candidateIsDominated                = all(Y <= C,1) & diffTest;
            %
            if ~any(candidateIsDominated)
               F(end+1)          = i;
            end
            %
            candidatesMask(i)                   = 0;
            candidatesMask(candidateDominates)  = 0;
         end
      end
      %% paretoSort
      function [F, Ys] = paretoSort(Y)
         len      = size(Y,2);
         F        = cell(len,1);
         Ys       = cell(len,1);
         i        = 1;
         unsorted = 1:len;
         while ~isempty(unsorted)
            F{i}     = unsorted(genetic.population.mgroup.getNonDominatedIndividuals(Y(:,unsorted)));
            Ys{i}    = Y(:,F{i});
            unsorted = setdiff(unsorted, F{i});
            i        = i+1;
         end
         F  = genetic.tools.removeEmptyCell(F);
         Ys = genetic.tools.removeEmptyCell(Ys);
      end
   end
end
