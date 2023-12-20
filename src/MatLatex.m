
%% carful to run this script variables from related scripts have to be run first

%load Data.mat
load Data2.mat
SNRMat = squeeze(SNRMat);
%%
for i = 1:2
    % get latex representation of tables
    table = squeeze(SNRMat(i,:,:));
    file = strcat("../report/Red",thesholding(i));
    file = strcat(file,".tex");
    columnLabels = num2cell(compose("%5.2e",p));
    rowLabels = num2cell(wavelets);
    
    matrix2latex(table,file,'format', '%.2f','rowLabels', rowLabels, 'columnLabels', columnLabels,'size', 'tiny');
end