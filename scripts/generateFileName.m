function filename = generateFileName(targetBox, c_dCas, plasmidCopyNum, dCasDissoCoe)
    fileTitle = "LineageData";
    c_dCas_str = convertCharsToStrings(sprintf('dCas9-%s', num2str(c_dCas)));
    target_str = convertCharsToStrings(sprintf('target-%s', targetBox));
    if plasmidCopyNum == 0 && dCasDissoCoe == 1
        filename = join([fileTitle, c_dCas_str, target_str], '_');
    elseif plasmidCopyNum > 0 && dCasDissoCoe == 1
        plasmid_str = convertCharsToStrings(sprintf('plasmid-%d', plasmidCopyNum));
        filename = join([fileTitle, plasmid_str, c_dCas_str, target_str], '_');
    elseif dCasDissoCoe ~= 1 && plasmidCopyNum == 0
        dCasDisso_str = convertCharsToStrings(sprintf('dCasOFF-%d', dCasDissoCoe));
        filename = join([fileTitle, c_dCas_str, target_str, dCasDisso_str], '_');
    else
        filename = 'NotExist';
    end
    filename = strcat(strrep(filename,'.','#'), '@*.mat');
    
    
    