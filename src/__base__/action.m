function output = action(list,varargin)

CommandSymbol = '*';
fprintf('* Please, make a selection: \n');

for i = 1:length(list)
   if isempty(varargin)
       cprintf('Keyword',['\t','  ','[',num2str(i),'] ',list{i}]);
   else
       cprintf('Keyword',['\t','  ',list{i}]);
   end
   fprintf('\n'); 
end

fprintf(CommandSymbol);
fprintf('  '); 

cprintf('Text','>> ');
if isempty(varargin), output = input('');
else output = input('','s');
end
end


