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
function methods()
% METHODS
%
%  Mono objective
%
%     Evolutionary/Stochastic
%     <a href="matlab:help genetic.optimizer.mo.annealing">annealing</a>   - Simulated Annealing algorithm
%     <a href="matlab:help genetic.optimizer.mo.cmaes">cmaes</a>       - Covariance Matrix Adaptation Evolution Strategy 
%     <a href="matlab:help genetic.optimizer.mo.pso">pso</a>         - Particle Swarm Optimisation algorithm
%
%     Direct Search Methods
%     <a href="matlab:help genetic.optimizer.mo.gss">gss</a>         - Generating Set Search method 
%     <a href="matlab:help genetic.optimizer.mo.mdsearch">mdsearch</a>    - Multi-Directional Search 
%     <a href="matlab:help genetic.optimizer.mo.simplex">simplex</a>     - Nelder-Mead Simplex algorithm
%
%     Local descent methods (use gradient)
%     <a href="matlab:help genetic.optimizer.mo.linesearch">linesearch</a>  - Line search method for local optimisation
%
%  Multi objective
%     <a href="matlab:help genetic.optimizer.mu.nsga2">nsga2</a>       - Non-dominated Sorting Genetic Algorithm II 
%     <a href="matlab:help genetic.optimizer.mu.nsga3">nsga3</a>       - Non-dominated Sorting Genetic Algorithm III
%     <a href="matlab:help genetic.optimizer.mu.spea2">spea2</a>       - Strength Pareto Evolutionary Algorithm II
%     <a href="matlab:help genetic.optimizer.mu.pesa2">pesa2</a>       - Pareto Envelope-based Selection Algorithm II
%     <a href="matlab:help genetic.optimizer.mu.mopso">mopso</a>       - Multi-Objective Particle Swarm Optimization algorithm 
%     <a href="matlab:help genetic.optimizer.mu.mombi2">mombi2</a>      - Many-Objective Metaheuristic Based on the R2 Indicator II
%     <a href="matlab:help genetic.optimizer.mu.moeadd">moeadd</a>      - Multi-Objective Evolutionary Algorithm based on Dominance and Decompositon
%     <a href="matlab:help genetic.optimizer.mu.tdea">tdea</a>        - Theta-Dominance based Evolutionary Algorithm
%     <a href="matlab:help genetic.optimizer.mu.rveastar">rveastar</a>    - Reference Vector (RV) guided Evolutionary Algorithm with RV regeneration strategy
help genetic.methods
end