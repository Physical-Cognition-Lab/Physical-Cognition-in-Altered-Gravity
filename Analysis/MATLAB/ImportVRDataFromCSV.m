clear; clc; close all;
vr_file_data = csvimport('FlexGravityVR.csv');
vr_data = vr_file_data(2:end,[2 11 12 15 16 18:21]);

% convert planet column
vr_data(cell2mat(cellfun(@(elem) strcmp(elem, 'b'), vr_data, 'UniformOutput', false))) = {-1};
vr_data(cell2mat(cellfun(@(elem) strcmp(elem, 'd'), vr_data, 'UniformOutput', false))) = {-1};
vr_data(cell2mat(cellfun(@(elem) strcmp(elem, 'w'), vr_data, 'UniformOutput', false))) = {-1};

vr_data(cell2mat(cellfun(@(elem) strcmp(elem, 'n'), vr_data, 'UniformOutput', false))) = {0};
vr_data(cell2mat(cellfun(@(elem) strcmp(elem, 'e'), vr_data, 'UniformOutput', false))) = {1};
vr_data(cell2mat(cellfun(@(elem) strcmp(elem, 'm'), vr_data, 'UniformOutput', false))) = {2};
vr_data(cell2mat(cellfun(@(elem) strcmp(elem, 'j'), vr_data, 'UniformOutput', false))) = {3};
vr_data(cell2mat(cellfun(@(elem) strcmp(elem, 'v'), vr_data, 'UniformOutput', false))) = {4};
vr_data(cell2mat(cellfun(@(elem) strcmp(elem, 'p'), vr_data, 'UniformOutput', false))) = {5};


% convert successful column
vr_data(cell2mat(cellfun(@(elem) strcmp(elem, 'y'), vr_data, 'UniformOutput', false))) = {1};

% convert distance column
vr_data(cell2mat(cellfun(@(elem) strcmp(elem, '.'), vr_data, 'UniformOutput', false))) = {-1};

% convert hand column
vr_data(cell2mat(cellfun(@(elem) strcmp(elem, 'l'), vr_data, 'UniformOutput', false))) = {1};
vr_data(cell2mat(cellfun(@(elem) strcmp(elem, 'r'), vr_data, 'UniformOutput', false))) = {0};

vr_data(cell2mat(cellfun(@(elem) isempty(elem), vr_data, 'UniformOutput', false))) = {-10};


vr_data = cellfun(@(elem) GetElemNumber(elem) , vr_data);

save('FlexGravityVR.mat','vr_data');




function elem_out=GetElemNumber(elem)
if (~isnumeric(elem))
    elem_out=str2num(elem);
else
    elem_out=elem;
end

end