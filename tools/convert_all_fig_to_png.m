function [  ] = convert_all_fig_to_png( pathN )
%CONVERT_ALL_FIG_TO_PNG will convert all the fig files in a folder to png
%files

if nargin > 0
    fig_files=dir(fullfile(pathN,'*.fig'));
    for i=1:length(fig_files)
    open(fullfile(pathN,fig_files(i).name));
    print('-dpng', '-r600', strcat(fig_files(i).name(1:end-4),'.png'));
    close all;
    end
else
    % Searching using wildcards
    fig_files=dir(fullfile(pwd,'*.fig'));
    for i=1:length(fig_files)
    open(fig_files(i).name);
    print('-dpng', '-r600', strcat(fig_files(i).name(1:end-4),'.png'));
    close all;
    end
end 

 
end