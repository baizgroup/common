function FixFTIRFileName (Path, SampleName)
%   This function renames the FTIR data files           %
%	For example: file '*.1' will be renamed to '*.01'   %
%   An already existing file will NOT be overwritten    %
%   Converted files will return an error message        %
%   Baiz Group. Fall 2019                               %
%   2019.10.7	Sherry.Y.                               %
%   Updated 2019.10.                                    %

%   Direct to the data folder
%   Path='D:\Baiz Group\Data\FTIR\20191003_MeSCN_temp\'

for k=0:9
    %  file1=[folder sprintf('1M_MeSCN_400_Dectrose-D-Glu.%d ',k)];
    %  file2=[folder sprintf('1M_MeSCN_400_Dectrose-D-Glu.0%d',k)]
    file1=[Path SampleName '.' sprintf('%d ',k)];
    file2=[Path SampleName '.' sprintf('0%d',k)];
    if isfile(file1) & ~isfile(file2)
        % File exists.
        movefile(file1 ,file2);
    elseif ~isfile(file1)
        % File1 does not exist.
        fprintf('Warning: file does not exist:\n%s', file1 );
    elseif isfile(file2)
        % File2 already exists.
        fprintf('Warning: file exists, can not overwrite:\n%s', file2 );
    end
 end