function curves = plotTable(A, col1, col2, f)
    if nargin < 4, f = @(x) x; end
    xData = eval(strcat(getVarStr(A),'.',col1)); 
    newcol2 = cell(size(col2));
    curves = cell(size(col2));
    for i = 1:length(col2)
        newcol2{i} = strrep(col2{i},'_',' ');
        yData = eval(strcat(getVarStr(A),'.',col2{i}));
        curves{i} = plot(xData(xData~=0),f(yData(xData~=0)));
        hold on
    end
%     l = legend(newcol2);
%     hold off
%     set(l, "Location", "Best")
    function varstr = getVarStr(~)
        varstr = sprintf('%s',inputname(1));
    end
end