function unpack_struct(s)
% Unpack structure by creating variables in the caller workspace 
% with names specified by field names
Field_Names = fieldnames(s);
for i = 1:length(Field_Names)
    assignin('caller',Field_Names{i}, s.(Field_Names{i}));
end