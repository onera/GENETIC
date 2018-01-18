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
function tab = tabular(T, headers, entries, options)
% Returns a formated tabular.
%
% Calling Sequence
%    tab = tabular(title, headers, entries)
%    tab = tabular(title, headers, entries, options)
%
% Parameters
%    title    : title of the table (string)
%    headers  : headers of each column (cell or matrix)
%    entries  : entries of the table (cell or matrix)
%    options  : options for the formating (structure)
%    tab      : formated tabular (string)
%
% Description
%    This routine enables to create formated tables.
%
% Example
%  A        = rand(10,3);
%  data     = [mean(A,2),max(A,[],2),min(A,[],2)];
%  T        = 'Table of results';
%  headers  = {'Mean','Max','Min'};
%  opt.float_format  = '%.3f';
%  opt.count_lines   = false;
%  opt.table_format  = 'simple';
%  tab               = tabular(T, headers, data, opt);
%  sprintf(tab)
%

% Options handling .........................................................
if nargin < 4
   options = struct();
end
default_options   = struct('float_format' ,  '%.2f',...
                           'table_format' ,  'simple',...
                           'alignment'    ,     [],...
                           'padding'      ,     1,...
                           'count_lines'  ,     false);
%
default_fields    = fieldnames(default_options);
opt_fields        = fieldnames(options);
tabular_check_for_wrong_opt(opt_fields, default_fields);
options           = tabular_complete_options(options, default_options);
% Reformating inputs .......................................................
% Entries of the tab
entries = tabular_format_entries(entries);
if options.count_lines
   for i = 1:length(entries)
      entries{i} = [sprintf('%d',i), entries{i}];
   end
end
entries = tabular_stringify(entries, options.float_format);
% Headers of the tab
headers = tabular_format_headers(headers);
if options.count_lines && ~isempty(headers)
   headers{1} = [' ', headers{1}];
end
headers = tabular_stringify(headers, options.float_format);
% Counting stuff ...........................................................
ncols       = length(entries{1});
if isempty(options.alignment)
   options.alignment = {};
   for i=1:ncols
      options.alignment{end+1} = 'c';
   end
   if options.count_lines
      options.alignment{1} = 'r';
   end
else
   if options.count_lines
      options.alignment = ['r',options.alignment];
   end
end
symbols.npadding     = options.padding;
[tab_len, cols_len]  = tabular_count_len(T, [headers,entries], ncols, options.padding);
% Preparing format .........................................................
symbols.lver         = '';
symbols.rver         = '';
symbols.mver         = '';
symbols.rdec         = '';
symbols.ldec         = '';
tabular_dec_text     = @(x) false;
symbols.dec          = tabular_dec_text;
top_bar              = '';
in_bar               = '';
bot_bar              = '';
switch options.table_format
   case 'simple'
      symbols.corner       = '+';
      symbols.lver         = '|';
      symbols.rver         = '|';
      symbols.mver         = '|';
      symbols.hor          = '-';
      top_bar              = tabular_make_line_sep(tab_len, symbols);
      bot_bar              = top_bar;
      head_bar             = tabular_make_line_sep(tab_len, symbols, cols_len);
      in_bar               = head_bar;
   case 'github-md'
      if ~isempty(T)
         disp(sprintf('''%s'' formating does not accept a title. Discarding.\n', options.table_format));
         T = [];
      end
      symbols.corner       = '|';
      symbols.lver         = '|';
      symbols.rver         = '|';
      symbols.mver         = '|';
      symbols.hor          = '-';
      head_bar             = tabular_make_line_sep(tab_len, symbols, cols_len, options.alignment);
   case 'latex'
      al = options.alignment;
      if options.count_lines
         al = [al{1},{'|'}, al{2:end}];
      end
      top_bar              = ['\\begin{table}\n\\begin{tabular}{',strcat(al{:}), '}\n'];
      head_bar             = '\\hline\n';
      bot_bar              = '\\end{tabular}\n';
      symbols.rver         = '\\\\';
      symbols.mver         = '&';
      symbols.rdec         = '$';
      symbols.ldec         = '$';
      symbols.dec          = @(x)tabular_isnum(x);
      if ~isempty(T)
         bot_bar = [bot_bar, '\\caption{', T, '}\n'];
      end
      bot_bar = [bot_bar, '\\end{table}'];
      T = '';
   otherwise
      error('Unknown table format');
end
% Building the tabular .....................................................
tab = top_bar;
if ~isempty(T)
   title_len   = tab_len - 2 - 2 * symbols.npadding;
   tab         = [tab, tabular_make_line({T}, [title_len], {'c'}, symbols)];
   tab         = [tab, top_bar];
end
if ~isempty(headers)
   tab = [tab, tabular_make_line(headers{1}, cols_len, options.alignment, symbols)];
   tab = [tab, head_bar];
end

for i = 1:length(entries)
   line = entries{i};
   if ~isempty(line)
      tab = [tab, tabular_make_line(line, cols_len, options.alignment, symbols)];
   else
      tab = [tab, in_bar];
   end
end
tab = [tab, bot_bar];
end

function out = tabular_justify(text, side, n)
text     = deblank(text);
n_text   = length(text);
d        = n - n_text;
switch side
   case 'l'
      out   = [text, tabular_blanks(d)];
   case 'r'
      out   = [tabular_blanks(d), text];
   case 'c'
      if mod(d,2) == 0
         n_l = d/2;
         n_r = d/2;
      else
         n_l = floor(d/2);
         n_r = ceil(d/2);
      end
      out = [tabular_blanks(n_l), text, tabular_blanks(n_r)];
end
end

function out = tabular_blanks(n)
blank = ' ';
out   = tabular_rep(blank, n);
end

function out = tabular_rep(symbol, n)
out = '';
for i = 1:n
   out = [out, symbol];
end
end

function [tab_len, cols_len] = tabular_count_len(T, entries, ncols, npadding)
if ~isempty(T)
   title_len = length(deblank(T));
else
   title_len = 0;
end
cols_len = zeros(ncols,1);
for i = 1:length(entries)
   line = entries{i};
   for j = 1:length(line)
      cols_len(j) = max(cols_len(j), length(line{j}));
   end
end
tab_len  = sum(cols_len) + 2*npadding * ncols + 1 * (ncols+1);
N        = title_len + 2  + 2;
if N > tab_len
   rem      = mod(N - tab_len, ncols);
   cols_len = cols_len + floor((N - tab_len) / ncols);
   tab_len  = title_len + 2 + 2 * npadding;
   if rem ~= 0
      cols_len = cols_len + 1;
      tab_len  = tab_len + (ncols-rem);
   end
end
end

function line = tabular_make_line(entries, cols_len, justif, symbols)
line = symbols.lver;
for i = 1:length(entries)
   text  = tabular_justify(entries{i}, justif{i}, cols_len(i));
   if symbols.dec(text)
      text = [symbols.ldec, text, symbols.rdec];
   end
   line = [line, tabular_add_padding(text, symbols.npadding)];
   if i < length(entries)
      line = [line, symbols.mver];
   end
end
line = [line, symbols.rver, '\n'];
end

function str = tabular_add_padding(str, n)
blank = ' ';
str   = [tabular_rep(blank,n), str, tabular_rep(blank, n)];
end

function line = tabular_make_line_sep(line_len, symbols, cols_len, alignment)
if nargin == 2
   n_corner = 2;
   line     = [symbols.corner, tabular_rep(symbols.hor, line_len - n_corner), symbols.corner];
elseif nargin == 3
   line = symbols.corner;
   for i = 1:length(cols_len)
      line = [line, tabular_rep(symbols.hor, cols_len(i) + 2*symbols.npadding)];
      if i < length(cols_len)
         line = [line, symbols.corner];
      end
   end
   line = [line, symbols.corner];
elseif nargin == 4
   line = symbols.corner;
   for i = 1:length(cols_len)
      switch alignment{i}
         case 'l'
            pre   = ':';
            post  = '';
         case 'c'
            pre   = ':';
            post  = pre;
         case 'r'
            pre   = '';
            post  = ':';
         otherwise
            pre   = '';
            post  = '';
      end
      add_size = length(pre) + length(post);
      line     = [line, pre, tabular_rep(symbols.hor, cols_len(i) + 2*symbols.npadding -add_size), post];
      if i < length(cols_len)
         line = [line, symbols.corner];
      end
   end
   line = [line, symbols.corner];
end
line = [line, '\n'];

end

function tab = tabular_stringify(tab, float_format)
if isempty(tab)
   return
end
for i=1:length(tab)
   line = tab{i};
   for j = 1:length(line)
      e = line{j};
      if isa(e,'double')
         tab{i}{j} = sprintf(float_format, e);
      end
   end
end
end

function tab = tabular_mat_to_tab(mat)
tab = {};
for i = 1:size(mat,1)
   inner_list = {};
   for j = 1:size(mat,2)
      inner_list{end+1} = mat(i,j);
   end
   tab{i} = inner_list;
end
end

function entries = tabular_format_entries(entries)
if isa(entries, 'double')
   % If the input is a matrix, then it is turned into a list  which contains
   % each row as a list
   entries  = tabular_mat_to_tab(entries);
elseif isa(entries, 'cell')
   % if the input is already a list, then one goes through all its elements
   % (meant to represent lines) to turn matrix into list of elements.
   for i = 1:length(entries)
      line = entries{i};
      if isa(line, 'double')
         inner_list = {};
         for j = 1:size(line, 2)
            inner_list{end+1} = line(j);
         end
         entries{i} = inner_list;
      end
   end
end
end

function headers = tabular_format_headers(headers)
if isempty(headers)
   return
end
if isa(headers, 'double')
   headers = tabular_mat_to_tab(headers);
elseif isa(headers, 'cell')
   headers = {headers};
end
end

function tabular_check_for_wrong_opt(fields, def_fields)
wrong_opt   = setdiff(fields, def_fields);
if ~isempty(wrong_opt)
   warn_msg       = strcat(wrong_opt',''', ''');
   disp('Warning: option(s) '''+warn_msg+''' are invalid. Discarding.');
   possible_opt   = strcat(def_fields',''', ''');
   disp('Possible values are: '''+possible_opt+'''.');
end
end

function options = tabular_complete_options(options, default_options)
default_fields    = fieldnames(default_options);
for f = default_fields'
   f = f{1};
   if ~isfield(options, f)
      options.(f) =  default_options.(f);
   end
end
end

function r = tabular_isnum(a)
if ( isnumeric(a) )
   r = 1;
else
   o = str2double(a);
   r = ~isempty(o);
end
end
 