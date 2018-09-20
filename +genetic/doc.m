function doc(page)
% GENETIC.DOC - Open the html documentation of the toolbox
%
% Syntax
%  genetic.doc

if nargin ==0
    page = 'index.html';
end
    p = genetic.tools.getPath();
    p = strrep(p,'+genetic','html');
    p = [p,filesep(),page];
    web(p);
end