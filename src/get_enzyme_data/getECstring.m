function EC_set = getECstring(EC_set,ecNumbers)
% getECstring  Build a single string of EC numbers for a protein.
%
% Provides a single string with all the EC numbers associated to a given
% protein, according to the following format:
% "ECX.X.X.X ECX.X.X.X ECX.X.X.X".
%
% Parameters
% ----------
% EC_set : char
%     existing EC string to which the EC numbers are appended.
% ecNumbers : char
%     EC numbers associated to the protein, separated by white space and/or
%     semicolons.
%
% Returns
% -------
% EC_set : char
%     a single string with all the EC numbers, each prefixed with "EC".
%
% Examples
% --------
%     EC_set = getECstring(EC_set, ecNumbers);
%
% See also
% --------
% findECInDB
%In case of several ec numbers for the same protein
new_EC_set = regexp(strtrim(ecNumbers),'[;\s]+','split');
new_EC_set = new_EC_set(~cellfun(@isempty,new_EC_set));
for l = 1:length(new_EC_set)
    EC_set = [EC_set 'EC' new_EC_set{l}];
    if l<length(new_EC_set)
        EC_set = [EC_set ' '];
    end
end
end
