% GENETIC.BENCH.EXPERIMENT is the interface for benchmarking methods.
%
% Syntax
%  xp = genetic.bench.experiment({{benchName_1, benchDim_1},...{benchName_N,  benchDim_N}},...
%                                 {method_1, ...., method_M})
%  xp = genetic.bench.experiment({{benchName_1, benchDim_1}, ..., {benchName_N,  benchDim_N}},...
%                                {{method_1, opt_1}        , ..., {method_M, opt_M}})
% Parameters
%  benchName_i : name of the benchmark
%  benchDim_i  : x-dimension of the benchmark
%  method_i    : name of the optimisation method
%  opt_i       : option for the corresponding optimisation method
%
% Description
%  This class is used to define experiments aimed a benchmarking
%  optimisation methods. More specifically, an experiment gathers a set
%  of benchmark functions, a set of optimisation methods (and associated
%  options).
%  Launching the experiment will apply each method on each benchmark,
%  store the results and compute some statistical indicators to ease the
%  interpretation of the raw data.
%
% Example
%  benchs   = {{'ackley',2}, {'griewank',2}};
%  methods  = {'cmaes','pso'};
%  xp       = genetic.bench.experiment(benchs, methods);
%  xp.name  = 'example';
%  xp.start();
%

% -------------------------------------------------------------------------

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
classdef experiment < handle
   properties
      benchList   = {};
      methodList  = {};
      optList     = {};
      %
      nRuns       = 10;
      %
      results     = [];
      target      = '';
      name        = '';
      stats       = [];
   end % end of properties
   methods
      function self = experiment(bench, methods, opt)
         if isa(bench, 'char') || isa(bench,'cell') && isa(bench{1}, 'char') && isa(bench{2}, 'double')
            bench = {bench};
         end
         for i = 1:length(bench)
            b = bench{i};
            if isa(b,'char')
               self.benchList{end+1}      = genetic.bench.load(b);
            elseif isa(b, 'cell') && isa(b{1}, 'char') && isa(b{2}, 'double')
               self.benchList{end +  1}   = genetic.bench.load(b{1}, b{2});
            else
               error('Expecting strings (and integers) for the bench.')
            end
         end
         %
         if isa(methods,'char')  || isa(methods, 'cell') && isa(methods{1},'char') && isa(methods{2},'struct')
            methods = {methods};
         end
         for i = 1:length(methods)
            m = methods{i};
            if isa(m, 'char')
               self.methodList{end+1}  = m;
               self.optList{end+1}     = [];
            elseif isa(m, 'cell') && isa(m{1},'char') && isa(m{2},'struct')
               self.methodList{end+1}  = m{1};
               self.optList{end+1}     = m{2};
            else
               error('Expecting strings (and structs) for the methods.');
            end
         end
         if nargin == 3 && ~isempty(opt)
            self.fillAttributes(opt);
         end
         if isempty(self.name)
            self.name = strrep(datestr(now),' ','-');
         end
      end
      %
      function fillAttributes(self, options)
         if ~isempty(options)
            fieldsNames = fieldnames(options);
            for i = 1:length(fieldsNames)
               try
                  self.(fieldsNames{i}) = options.(fieldsNames{i});
               catch
                  fprintf('''%s'' is not an attribute of %s, discarding.\n',fieldsNames{i},class(self));
               end
            end
         end
      end
      %
      function start(self, force)
         % GENETIC.BENCH.EXPERIMENT.START starts the experiment, computes
         % statistical data with the results and stores everything.
         if nargin < 2
            force = false;
         end
         if isempty(self.target) || force
            % target file
            self.target = [genetic.bench.experiment.getPath(),self.name];
         end
         if ~force && exist([self.target,'.mat'], 'file')
            fprintf('Target file already exists, use force option to replace it.\n');
            return
         end
         self.results   = struct();
         ntot           = length(self.benchList) * length(self.methodList) * self.nRuns;
         iter           = 0;
         for i = 1:length(self.benchList)
            % for each bench ...
            b           = self.benchList{i};
            [fun,n,cst] = genetic.bench.toCall(b);
            for j = 1:length(self.methodList)
               % for each method ...
               method   = self.methodList{j};
               opt      = self.optList{j};
               R        = cell(self.nRuns,1);
               for k = 1:self.nRuns
                  % for each run ...
                  iter  = iter+1;
                  fprintf('Performing experiment...%.1f%%\n',iter/ntot*100);
                  res      = genetic.bench.experiment.resultStructure();
                  res.xDim = b.xDim;
                  res.fDim = b.fDim;
                  try
                     % one tries to launch the optimisation and to store
                     % the results (as well as various indicators)
                     [xopt, fopt, info]   = genetic.min(fun, n, method, cst, opt);
                     res.success          = info.success;
                     res.error            = false;
                     res.xopt             = xopt;
                     res.fopt             = fopt;
                     res.flag             = info.stopFlags;
                     res.nImprove         = info.nImprove;
                     res.nEval            = info.nEval;
                     res.elapsedTime      = info.elapsedTime;
%                      res.maxCstViolation  = info.maxCstViolation;
                  catch ME
                     % if the optimisation fails (with an error), then one
                     % stores the error.
                     res.success          = false;
                     res.error            = true;
                     res.errorMsg         = ME;
                  end
                  R{k} = res;
               end
               % The results are stored in a structure so that:
               %     self.results.benchName.methodName
               % gathers the results (one for each run)
               self.results.(b.name).(method) = R;
            end
         end
         % Once all the optimisation have been performed, some statistical
         % indicators are computed on the results
         self.makeStats();
         % Then, the object itself is saved.
         save(self.target,'self');
         fprintf('Experiment saved at ''%s''\n',self.target);
      end
      %
      function makeStats(self)
         if isempty(self.results)
            fprintf('Statistics can only be created when the bench have been launched.\n');
            return
         end
         R         = self.results;
         for i = 1:length(self.benchList)
            bench = self.benchList{i};
            for j = 1:length(self.methodList)
               method   = self.methodList{j};
               rij      = R.(bench.name).(method);
               if bench.fDim == 1
                  sta      = genetic.bench.stats(rij);
               else
                  error('Multi-objective statistics not yet implemented.\n');
               end
              
               self.stats.(bench.name).(method) = sta;
            end
         end
      end
      %
      function showTab(self, tabType, refs, field)
         if isa(refs,'char')
            refs = {refs};
         end
         %
         S           = self.stats;
         switch tabType
            case 'm'
               benchNames  = fieldnames(S);
               for i = 1:length(refs)
                  method   = refs{i};
                  entries  = cell(0);
                  for b = benchNames'
                     b              = b{1};
                     res            = S.(b).(method).(field);
                     bName          = [b,' (',num2str(self.results.(b).(method){1}.xDim),')'];
                     entries{end+1} = {bName, res.min, res.mean, res.max, res.std};
                  end
                  title    = [method,' (',field,')'];
                  headers  = {'bench','min','mean','max','std'};
                  opt.float_format = '%1.2e';
                  tab      = genetic.tools.tabular(title, headers, entries, opt);
                  fprintf(tab);
               end
            case 'b'
               methodNames  = fieldnames(S.(refs{1}));
               for i = 1:length(refs)
                  b        = refs{i};
                  entries  = cell(0);
                  for m = methodNames'
                     m              = m{1};
                     res            = S.(b).(m).(field);
                     mName          = m;
                     entries{end+1} = {mName, res.min, res.mean, res.max, res.std};
                  end
                  title    = [b,' (',num2str(self.results.(b).(m){1}.xDim),' - ',field,')'];
                  headers  = {'method','min','mean','max','std'};
                  opt.float_format = '%1.2e';
                  tab      = genetic.tools.tabular(title, headers, entries, opt);
                  fprintf(tab);
               end 
         end
      end
      %
   end % end of methods
   methods (Static)
      function out = resultStructure()
         out = struct('xDim'     , [],...
                      'fDim'     , [],...
                      'error'    , [],...
                      'errorMsg' , [],...
                      'success'  , [],...
                      'xopt'     , [],...
                      'fopt'     , [],...
                      'flag'     , [],...
                      'nEval'    , [],...
                      'nImprove' , [],...
                      'elapsedTime',[],...
                      'maxCstViolation',[]);
      %
      end
      %
      function xpPath = getPath()
         info        = what('+genetic');
         geneticPath = info.path;
         xpPath      = [geneticPath,filesep,'+bench',filesep,'xp',filesep];
         if ~isdir(xpPath)
            mkdir(xpPath);
         end
      end
      %
      function xp = load(name)
         xpPath   = genetic.bench.experiment.getPath();
         load([xpPath,name,'.mat']);
         xp       = self;
      end
   end
end