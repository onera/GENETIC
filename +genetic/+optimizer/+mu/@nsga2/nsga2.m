classdef nsga2 < genetic.optimizer.multi & genetic.optimizer.simpleScheme
   % NSGA2 (Non-dominated Sorting Genetic Algorithm II)
   %
   %  NSGA2 is by far the most well known and most used multi-objective
   %  optimization metaheuristic.
   %  
   %  NSGA2 is an efficient multi-objective evolutionary algorithm based on
   %  an elitist approach. Its particular fitness assignment consists in
   %  sorting the population in different fronts using the non-domination
   %  order relation. To form the next generation, the algorithm combines
   %  both current population and corresponding offspring, generated with
   %  the standard differential evolution-based or bimodal crossover or
   %  polynomial-based mutation operators. Finally, the best individuals in
   %  terms of non-dominance and diversity are chosen.
   %
   %  NSGA2 has a low time complexity of O(N x log(N)), where integer N is
   %  the population size.
   %
   % References
   %  [1] K. Deb, Multi-Objective Optimization Using Evolutionary Algorithms.
   %      John Wiley & Sons, LTD, 2001.
   %  [2] K. Deb, A. Pratap, S. Agarwal and T. Meyarivan, A fast and elitist
   %      multiobjective genetic algorithm: NSGA-II. In IEEE Transactions on
   %      Evolutionary Computation, vol. 6, no. 2, pp. 182-197, Apr 2002.
   %  [3] A.J. Nebro, J.J. Durillo, M. Machin, C.A. Coello Coello, B. Dorronsoro,
   %      A Study of the Combination of Variation Operators in the NSGA-II Algorithm.
   %      Proceedings of the 15th Conference of the Spanish Association for Artificial
   %      Intelligence, CAEPIA 2013, Madrid, Spain, September 17-20, 2013. Lecture
   %      Notes in Computer Science Volume 8109, pp. 269-278, 2013.
   %
   properties
      sbxCrossoverProbability = 2/3;
      polyMutationProbability = 1/3;
      crossoverDistribIndex   = 20;
      mutationDistribIndex    = 20;
      CR                      = 0.8803;
      F                       = 0.4717;
      localGroup
   end
   properties
      X = [];
      Y = [];
   end
   methods
      %% Constructor
      function self = nsga2(varargin)
         self@genetic.optimizer.multi(varargin{:})
         self@genetic.optimizer.simpleScheme()
         self.methodName         = 'nsga2';
         self.longName           = 'Non-dominated Sorting Genetic Algorithm II';
         self.needFiniteBounds   = true;
      end
      %% Evolve
      function group = evolve(self, group)
         %
         if ~isempty(self.X)
            % Creating virtual group gathering the freshly evaluated points
            % and the ones that have been kept previously
            Yaug                    = [self.Y, group.getObjective()];
            Xaug                    = [self.X, group.getValue()];
            [rankedFronts, Yaugs]   = genetic.population.mgroup.paretoSort(Yaug);
            %
            nRanks                  = length(rankedFronts);
            % 
            CD = cell(nRanks,1);
            for i = 1:nRanks
               CD{i} = genetic.optimizer.mu.nsga2.crowdingDist(Yaugs{i});
            end
            % Selecting the best individuals in the virtual group
            idx                     = genetic.optimizer.mu.nsga2.selectBest(rankedFronts, nRanks, CD, group.len);
            % Updating P
            self.X                  = Xaug(:,idx);
            self.Y                  = Yaug(:,idx);
            %
            currentX                = self.X;
         else
            % First iteration
            self.X                  = group.getValue();
            self.Y                  = group.getObjective();
            currentX                = self.X;
         end
         % New values for the group
         newX              = self.crossover(currentX, group.len);
         group.moveTo(newX);
      end
      %% Crossover
      function y = crossover(obj, x, popSize)
         y = [];
         while size(y,2) < popSize
            randValue = rand(1);
            if randValue <= obj.polyMutationProbability
               individuals = genetic.population.group.select(1, x, 'rand');
               tmp         = genetic.population.group.polyMutation(individuals, obj.xMin, obj.xMax, obj.polyMutationProbability, obj.mutationDistribIndex);
            elseif randValue <= obj.sbxCrossoverProbability
               if popSize >= 2
                  individuals = genetic.population.group.select(2, x, 'rand');
                  tmp         = genetic.population.group.sbxCrossover(individuals, obj.xMin, obj.xMax, [], obj.crossoverDistribIndex);
               else
                  individuals = genetic.population.group.select(1, x, 'rand');
                  tmp         = genetic.population.group.polyMutation(individuals, obj.xMin, obj.xMax, obj.polyMutationProbability, obj.mutationDistribIndex);
               end
            else
               if popSize >= 3
                  [individuals, idxList]  = genetic.population.group.select(3, x, 'rand');
                  idx                     = randi(popSize);
                  while ~isempty(intersect(idxList,idx))
                     idx = randi(popSize);
                  end
                  candidateX  = x(:,idx);
                  tmp         = genetic.population.group.deCrossover(individuals, candidateX, obj.xMin, obj.xMax, obj.CR, obj.F);
               else
                  if popSize == 2
                     individuals = genetic.population.group.select(2, x, 'rand');
                     tmp         = genetic.population.group.sbxCrossover(individuals, obj.xMin, obj.xMax, [], obj.crossoverDistribIndex);
                  elseif popSize == 1
                     individuals = genetic.population.group.select(1, x, 'rand');
                     tmp         = genetic.population.group.polyMutation(individuals, obj.xMin, obj.xMax, obj.polyMutationProbability, obj.mutationDistribIndex);
                  end
               end
            end
            y = [y tmp];
            if size(y,2) > popSize
               y = y(:,1:popSize);
            end
         end
      end
   end
   %%
   methods (Static)
      function d = crowdingDist(Y)
         ny       = size(Y,1);
         n        = size(Y,2);
         [Ys,idx] = sort(Y, 2, 'ascend');
         %
         d        = zeros(ny,n);
         for i = 1:ny
            idi            = idx(i,:);
            d(i,idi(1))    = inf;
            d(i,idi(end))  = inf;
%             for j = 2:n-1
%                d(i,idi(j))      =  (Y(i,idi(j+1))-Y(i,idi(j-1)))/(Ys(i,end) - Ys(i,1));
%             end
            d(i,idi(2:end-1)) = (Y(i,idi(3:end)) - Y(i,idi(1:end-2))) ./ repmat(Ys(i,end) - Ys(i,1),1,n-2);
         end
         d = sum(d,1);
      end
      %% Select best points
      function idx = selectBest(rankedFronts, nRanks, CD, popSize)
         nX = 0;
         idx         = [];
         for i = 1:nRanks
            if (length(rankedFronts{i})+nX) <= popSize
               toAdd = rankedFronts{i};
            else
               n           = popSize - nX;
               [~,idxList] = sort(CD{i},'descend');
               toAdd       = rankedFronts{i};
               toAdd       = toAdd(idxList(1:n));
            end
            if ~isempty(toAdd)
               idx   = [idx toAdd];
            else
               break
            end
            nX = nX + length(toAdd);
         end
      end
   end
end