function [Data] = func_getStrainSelection(app)
    SelectedResults = get(app.StrainTree.SelectedNodes);
    if isempty(SelectedResults)
        Data = {}; % will throw error that no data is selected
    else
        % pair the selected images with its object
        for s = 1:length(SelectedResults)
            Child = SelectedResults(s).Text;
            Parent = SelectedResults(s).Parent.Text;

            % scan trhough the saved analysis list for matching Name
            for p = 1:size(app.SavedAnalysis,1)
                % if the parent node is the same as one of the listed
                % analysis objects return the object
                if isequal(app.SavedAnalysis{p,2},Parent)
                    objNum = p;
                end
            end

            Data{s,1} = str2num(Child(end-2:end));
            Data{s,2} = objNum;
            Data{s,3} = Child;
            Data{s,4} = Parent;

        end
    end

end

