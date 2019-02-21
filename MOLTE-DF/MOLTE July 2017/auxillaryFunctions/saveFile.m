function saveFile( name, choice, policies, problemi)
%A handful function for saving .mat file in a specific folder
   save(fullfile('./',int2str(problemi),name), 'choice','policies');  

end

