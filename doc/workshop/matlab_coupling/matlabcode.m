% test code for forward sampling of parameters with RAVEN

clear
clc
close all



var_test = $RAVEN-var_test$;
calculated_var = var_test*2;

% write csv
cHeader = {'var_test' 'calculated_var'}; %dummy header
commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas
% % write header to file
fid = fopen('matlab_output.csv','w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);
% write data to end of file
data = [var_test calculated_var];
dlmwrite('matlab_output.csv',data,'-append');



exit 
