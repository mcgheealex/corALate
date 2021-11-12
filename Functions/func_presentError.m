function [app] = func_presentError(app,error)
%FUNC_PRESENTERROR Summary of this function goes here
%   Detailed explanation goes here
    
    switch error
        case 'no strain data selected'
            disp('no strain data selected')
        otherwise
            disp('Code error and real error (•ω•`)')
    end
    
end

