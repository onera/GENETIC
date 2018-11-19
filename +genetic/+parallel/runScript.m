function runScript()
tmpPath              = genetic.tools.getPath('tmp');
data                 = load([tmpPath 'data.mat']);
[xopt, fopt, info]   = genetic.min(data.fun, data.x, data.method, data.constraints, data.options);
%
if ~strcmp(tmp.options.parallel.target,'local')
   exit;
end
end