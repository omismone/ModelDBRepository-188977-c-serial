%% be sure that Enoise and Inoise matrices are in the workspace before running the script
% to have them run the simulation at https://github.com/ModelDBRepository/188977

% save Enoise
fileID = fopen('bin/Enoise.bin','wb');
fwrite(fileID,Enoise','double');
fclose(fileID);

%save Inoise
fileID = fopen('bin/Inoise.bin','wb');
fwrite(fileID,Inoise','double');
fclose(fileID);
