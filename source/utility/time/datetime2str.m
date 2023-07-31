function s = datetime2str(s)
    if isstruct(s)
        fields = fieldnames(s);
        for i = 1:numel(fields)
            field = fields{i};
            value = s.(field);
            if isstruct(value)
                s.(field) = datetime2str(value);
            elseif isa(value, 'datetime')
                s.(field) = cellstr(datestr(value));
            end
        end
    end
end