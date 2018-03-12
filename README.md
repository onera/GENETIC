![genetic logo](./logo.png)

**Genetic** is a Matlab/Octave toolbox for optimisation which gathers several mono and multi-objective optimisation algorithms.

# Table of Contents

* [Installation](#installation)
* [Getting started](#getting-started)
   * [Main interface for optimisation](#main-interface-for-optimisation)
   * [Available methods](#available-methods)
      * [Mono-objective](#mono-objective)
      * [Multi-objective](#multi-objective)
   * [Options](#options)
   * [Constraints](#constraints)
   * [Examples](#examples)
      * [Mono-objective](#mono-objective-1)
      * [Multi-objective](#multi-objective-1)
* [Benchmarking tools](#benchmarking-tools)
   * [Benchmarks](#benchmarks)
   * [Comparing methods](#comparing-methods)
* [Contact](#contact)



# Installation

**Genetic** should work
* with Matlab >= R2014b,
* with Octave >= 4.0 (some functionalities might be missing though).

To install it, download the last release and add the directory to your path in Matlab/Octave.

# Getting started

## Main interface for optimisation

All the optimisation methods of **Genetic** are called through the unified interface ``genetic.min`` as follows:

~~~ Matlab
[xopt, fopt, info] = genetic.min(f, x, method)
[xopt, fopt, info] = genetic.min(f, x, method, options)
[xopt, fopt, info] = genetic.min(f, x, method, constraints, options)
~~~

where the input arguments are:
* ``f``: the objective function. For mono-objective optimisation, it can be only a handle function and for multi-objective optimisation, it must be a cell containing a handle function and the number of objectives, i.e. ``f = {f_handle, fDim}``.
* ``x``: the dimension of the problem or an initial point/population.
* ``method``: the name of the method to be used.
* ``options``: a structure containing the parameters for the optimisation process.
* ``constraints``: a structure containing a description of the constraints. (*not yet functional*)

and the outputs are:
* ``xopt``: the best solution(s) found.
* ``fopt``: the associated objective or Parefo front.
* ``info``: a structure containing informations about the optimisation process.

*Remark: Similarly, the main interface for maximisation is ``genetic.max``.*


## Available methods

**Genetic** gathers both mono and multi-objective optimisation methods listed below. Note that contrary to what its name may suggest, the toolbox does not contain only evolutionary methods.

To get additional help on the methods (and references), type ``genetic.methods``.

### Mono-objective

|   Key           |          Name                                       |
| ---------------:|:----------------------------------------------------|
| ``annealing``   | Simulated Annealing algorithm (WIP)                 |
| ``pso``         | Particle Swarm Optimisation algorithm               |
| ``cmaes``       | Covariance Matrix Adaptation Evolution Strategy     |
| ``simplex``     | Nelder-Mead Simplex algorithm                       |
| ``mdSearch``    | Multi-Directional direct Search                     |
|   ``gss``       | Generating Set Search Method                        |
| ``linesearch``  | Local descent algorithm based on line-search        |

### Multi-objective

|   Key           |          Name                                                                     |
| ---------------:|:----------------------------------------------------------------------------------|
| ``nsga2``       | Non-dominated Sorting Genetic Algorithm II                                        |
| ``nsga3``       | Non-dominated Sorting Genetic Algorithm III                                       |
| ``spea2``       | Strength Pareto Evolutionary Algorithm II                                         |
| ``pesa2``       | Pareto Envelope-based Selection Algorithm II                                      |
| ``mopso``       | Multi-Objective Particle Swarm Optimization algorithm                             |
| ``mombi2``      | Many-Objective Metaheuristic Based on the R2 Indicator II                         |
| ``moeadd``      | Multi-Objective Evolutionary Algorithm based on Dominance and Decompositon        |
| ``tdea``        | Theta-Dominance based Evolutionary Algorithm                                      |
| ``rveastar``    | Reference Vector (RV) guided Evolutionary Algorithm with RV regeneration strategy |

## Options

The ``options`` structure contains parameters that let you tune the optimisation process. In particular, it enables to modify both general and method-specific parameters. The general parameters, i.e., those shared by all the methods, are the following:

| Name                  |                Description                                            |  Default |
|----------------------:|:----------------------------------------------------------------------|:--------:|
| ``maxFunEval``        | Maximum number of function evaluation                                 | ``5000`` |
| ``verbosity``         | Verbosity level of the algorithm                                      | ``0``    |
| ``popSize``           | Dimension of the population (for population-based methods)            | ``round(10 + 2*sqrt(n))``      |
| ``gradObj``           | Availability of the gradient as second output of the objective function | ``false`` |

The general parameters shared by all the mono-objective methods are the following:

| Name                  |               Description                                             |  Default |
|----------------------:|:----------------------------------------------------------------------|:--------:|
| ``targetY``           | Target value of the objective function                                | ``-inf`` |
| ``tolY``              | Minimum decrease between two consecutive best solutions               | ``1e-8`` |
| ``ntolY``             | Number of time the minimum objective decrease must be met before stopping      | ``1`` |
| ``tolX``              | Minimum distance between two consecutive best solutions               | ``1e-8`` |
| ``ntolX``             | Number of time the minimum distance must be met before stopping       | ``1`` |


For method-specific parameters, refer to the corresponding help page.


## Constraints

*/!\ Constraints are not yet handled*

Similarly to options, constraints are added through a structure containing specific fields depending on the constraints:

  * bounds are characterized by the fields ``xMin`` and ``xMax`` such that: <b>x</b><sub>min</sub> &le; <b>x</b> &le; <b>x</b><sub>max</sub>,
  * linear inequality constraints are characterized by the fields ``A`` and ``b`` such that: A<b>x</b> &le; <b>b</b>,
  * linear equality constraints are characterized by ``Aeq`` and ``beq`` such that: A<sub>eq</sub><b>x</b>=<b>b</b><sub>eq</sub>,
  * non-linear inequality constraints are given as an anonymous function ``c`` such that: <b>c</b>(<b>x</b>) &le; 0,
  * non-linear equality constraints are given as an anonymous function ``ceq`` such that: <b>c</b><sub>eq</sub>(<b>x</b>) = 0.



## Examples

### Mono-objective

~~~ Matlab
xDim                = 10;
f                   = @(x) norm(x)^2;
method              = 'cmaes';
options             = struct('verbosity', 3, 'maxFunEval', 3000);
[xopt, fopt, info]  = genetic.min(f, xDim, method, options);
~~~

### Multi-objective

~~~ Matlab
xDim                = 10;
A                   = rand(xDim,xDim);
b                   = rand(xDim,1);
fun                 = @(x) [norm(x)^2;norm(A*x - b)^2];
fDim                = 2;
f                   = {fun, fDim};
method              = 'nsga2';
options             = struct('verbosity', 3, 'maxFunEval', 10000);
[xopt, fopt, info]  = genetic.min(f, xDim, method, options);
~~~

# Benchmarking tools

## Benchmarks

*Load and use benchmarks.* **Genetic** is shipped with a set of mono and multi-objective academic benchmark problems gathered from the literature. These problems that can be listed with ``genetic.bench.list()``.

The data associated with a problem can then be accessed with ``genetic.bench.load(key, n)`` where ``key`` is the name of the benchmark and ``n`` is the dimension of the problem. Note that some problems have a fixed dimension that is displayed when listing the problems.

The benchmarks can also be directly called from **Genetic** by replacing the objective function in the call with the name of the benchmark, i.e. ``genetic.min(benchName, x, method)``.

*Add your own benchmarks to genetic.* You can add your own benchmark to **Genetic** through ``genetic.bench.create`` as follows,

~~~ Matlab
% Creation of the benchmark
ben.f    = @(x) norm(x);
ben.xDim = 10;
ben.fopt = 0;
ben.xopt = zeros(10,1);
ben.name = 'norm10';
genetic.bench.create(ben);
clear ben
% Re-loading the benchmark
ben      = genetic.bench.load('norm10')
~~~

It can then be used as the other benchmarks. Note that it will fail during execution if the handle function in the benchmark requires some other function that is not in the Matlab path.

## Comparing methods

**Genetic** contains also some elements to help you compare the performances of its methods. It is done through the creation of a ``genetic.bench.experiment`` that let you run several methods (or the same with different options) on several benchmarks (or the same at different dimensions). For instance, to compare ``cmaes`` and ``pso`` on the ackley and griewank functions,

~~~ Matlab
% Creating the experiment
benchs   = {{'ackley',2}, {'griewank',2}};
methods  = {'cmaes','pso'};
xp       = genetic.bench.experiment(benchs, methods);
xp.name  = 'example';
% Run the optimisations
xp.start();
~~~

All the raw results of the optimisations are stored in ``xp.results``. To ease the comparison, statistical elements can be displayed with ``xp.showTab``, e.g. to show the distance with the optimal value:

~~~
>> xp.showTab('b', 'ackley', 'fopt')
+----------------------------------------------------+
|                 ackley (2 - fopt)                  |
+----------------------------------------------------+
| method |   min    |   mean   |   max    |   std    |
+--------+----------+----------+----------+----------+
| cmaes  | 9.05e-09 | 1.35e-07 | 6.46e-07 | 2.25e-07 |
|  pso   | 3.97e-09 | 8.33e-08 | 2.85e-07 | 8.55e-08 |
+----------------------------------------------------+
~~~

Other elements like ``nEval``, ``elapsedTime``, ``nImprove`` may also be displayed. Similarly, the statistics for one method among all the benchmarks can be displayed using ``xp.showTab('m', 'cmaes', 'fopt')``.


# Contact

Feel free to contact us at: genetic [at] onera.fr
