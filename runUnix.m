%% write .bat to invoke java 
% write .bat and call .bat
% arg1: fullpath of java class,
% arg2: call java or perl function with arguments required
function [s,w] = runUnix(arg1, arg2)
fid = fopen('a.bat','w');
fprintf(fid, '#!/bin/sh\n%s\n%s\n', ['cd ', arg1], arg2);
fclose(fid);
[s,w] = unix('sh a.bat');
    