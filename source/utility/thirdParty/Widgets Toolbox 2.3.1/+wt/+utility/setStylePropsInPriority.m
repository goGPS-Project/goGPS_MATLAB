function setStylePropsInPriority(comps, propNames, value)
% For each given component, set the specified value to the first identified
% property. Properties specified in prioritized order.

% Copyright 2023 The MathWorks Inc.


% Define arguments
arguments
    comps matlab.graphics.Graphics
    propNames (1,:) string
    value
end

% Filter any invalid components
comps(~isvalid(comps)) = [];

% Track components that have been set
isDone = false(size(comps));

% Loop on each property
for thisProp = propNames

    % Does the current property exist in each component?
    needsSet = ~isDone & isprop(comps, thisProp);

    % Set as needed
    if any(needsSet)
        set(comps(needsSet), thisProp, value);
        isDone(needsSet) = true;
    end

    % Return early if complete
    if all(isDone)
        return;
    end

end %for
