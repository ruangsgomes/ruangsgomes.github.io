function save_data (Output_MOGA,GA,j,algorithm)

      
% output file
if j == 1 && strcmp(algorithm, 'MODE')
    fileID = fopen('output_EPR.txt','w');
    fprintf(fileID,'%12s\r\n','---------- MODEGA-EPR Toolbox output file ----------');
    fprintf(fileID,'%12s\r\n','algorithm:',[ algorithm]);
    fprintf(fileID,'%12s %6.2f\r\n','EPR terms =',GA.m+1);
    fprintf(fileID,'%12s\r\n','Best population');
    fprintf(fileID,'%6.2f %6.2f %6.2f\n',Output_MOGA.bestpop');
    fprintf(fileID,'%12s\r\n','Best parameters');
    fprintf(fileID,'%6.6f \n',Output_MOGA.best_par');
    fprintf(fileID,'%12s\r\n','Statistics');
    fprintf(fileID,'%12s %6.2f\r\n','R^2 =',Output_MOGA.R2);
    fprintf(fileID,'%12s %6.2f\r\n','RMSE =',Output_MOGA.RMSE);
    fprintf(fileID,'%12s %6.2f\r\n','MAE =',Output_MOGA.MAE);
    fprintf(fileID,'%12s %6.2f\r\n','r =',Output_MOGA.r);
    fprintf(fileID,'%12s %6.2f\r\n','SSE =',Output_MOGA.sse_min);
    fprintf(fileID,'%12s %6.2f\r\n','E_rel =',Output_MOGA.E_rel);
    fclose(fileID);
    
else
    fileID = fopen('output_EPR.txt','a');
    fprintf(fileID,'%12s\r\n','--------------------');
    fprintf(fileID,'%12s\r\n','algorithm:',[ algorithm]);
    fprintf(fileID,'%12s %6.2f\r\n','EPR terms =',GA.m+1);
    fprintf(fileID,'%12s\r\n','Best population');
    fprintf(fileID,'%6.2f %6.2f %6.2f\n',Output_MOGA.bestpop');
    fprintf(fileID,'%12s\r\n','Best parameters');
    fprintf(fileID,'%6.6f \n',Output_MOGA.best_par');
    fprintf(fileID,'%12s\r\n','Statistics');
    fprintf(fileID,'%12s %6.2f\r\n','R^2 =',Output_MOGA.R2);
    fprintf(fileID,'%12s %6.2f\r\n','RMSE =',Output_MOGA.RMSE);
    fprintf(fileID,'%12s %6.2f\r\n','MAE =',Output_MOGA.MAE);
    fprintf(fileID,'%12s %6.2f\r\n','r =',Output_MOGA.r);
    fprintf(fileID,'%12s %6.2f\r\n','SSE =',Output_MOGA.sse_min);
    fprintf(fileID,'%12s %6.2f\r\n','E_rel =',Output_MOGA.E_rel);
    fclose(fileID);
       
end
