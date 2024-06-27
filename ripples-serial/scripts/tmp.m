%writematrix(num2str(vE','%.16f, '),'C:\Users\mclab\Desktop\simone\thesis\ripples-serial\ripples-serial\bin\tmp.txt');

fileID = fopen('C:\Users\mclab\Desktop\simone\thesis\ripples-serial\ripples-serial\bin\tmp.txt','w');
fprintf(fileID,"%.14f, ",vE');
fclose(fileID);
