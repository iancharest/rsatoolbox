% For some unknown reason this is necessary if using parfor
% Will save a variable called 'x'.
function parsave(savePath, x)  %#ok<INUSD>
    save('-v7.3', savePath, 'x')
end
