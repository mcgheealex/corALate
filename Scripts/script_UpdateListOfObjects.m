d = app.DisplacementTree.Children;
d.delete;
s = app.StrainTree.Children;
s.delete;

for i=1:size(app.SavedAnalysis,1)
    Items{i} = ['Name: ', app.SavedAnalysis{i,2} , '     Analysis type: ', app.SavedAnalysis{i,4}];
    
    obj = app.SavedAnalysis{i,1};
    
    % displacement selection box
    D = app.DisplacementTree;
    parentD{i} = uitreenode(D,'Text',app.SavedAnalysis{i,2});
    % strain selection box
    S = app.StrainTree;
    parentS{i} = uitreenode(S,'Text',app.SavedAnalysis{i,2});
    
    if isempty(obj.ResultDisp)
        childD = uitreenode(parentD{i},'Text','no images analyzed');
        childS = uitreenode(parentS{i},'Text','no images analyzed');
    else
        for j=1:length(obj.ResultDisp)
            childD = uitreenode(parentD{i},'Text',['IM_' sprintf('%03d',j)]);
            childS = uitreenode(parentS{i},'Text',['IM_' sprintf('%03d',j)]);
        end
    end
    
end

% export selection box
app.SelectAnalysisListBox_Export.Items = Items;