function [xopt, fopt, info] = execute(fun, x, method, constraints, options)
paraOpt     = options.parallel;
cpuTarget   = paraOpt.target;
nCores      = paraOpt.nCores;
%
if strcmp(cpuTarget,'local')
   [~, out] = system('hostname');
   cpuName  = deblank(out);
else
   cpuName  = cpuTarget;
end
%
tmpPath     = genetic.tools.getPath('tmp');
save([tmpPath, 'data'],'fun','x','method','constraints','options');
switch cpuTarget
   case 'local'
      %
      if options.verbosity > 0
         verbosity   = ' > mpiLogFile';
      else
         verbosity   = '';
      end
      %
      command        = ['mpirun -n ' num2str(nCores) ' matlab -nojvm -nodesktop -nodisplay -r genetic.parallel.runScript', verbosity];
      %
      system(command);
      %
      results              = load([tmpPath 'results.mat']);
      xopt                 = results.xopt;
      fopt                 = results.fopt;
      info                 = results.info;
      %
      delete([tmpPath '*.mat']);
end
end