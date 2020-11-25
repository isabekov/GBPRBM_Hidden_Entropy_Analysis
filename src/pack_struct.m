function pack_struct(sname, vars)
% Create a structure with the name specified by the string "sname"
% with field name listed in "vars" cell and the corresponding variable
% values.
evalin('caller', [sname ' = [];']);
for i = 1:length(vars)
    evalin('caller', [sname '.' vars{i} ' = ' vars{i} ';']);
end