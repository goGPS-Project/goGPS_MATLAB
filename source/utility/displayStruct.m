function displayStruct(s, indent)
    % Use the function
    %displayStruct(export_data)
    if nargin < 2
        indent = '  ';
    end
    fields = fieldnames(s);
    for i = 1:numel(fields)
        field = fields{i};
        value = s.(field);
        if isstruct(value)
            fprintf('%s%s:\n', indent, field);
            displayStruct(value, [indent '  ']);
        else
            fprintf('%s%s: %s, %s\n', indent, field, class(value), mat2str(size(value)));
        end
    end
end

