function writeInFile(file,N)
    fd= mopen(file,'a+');
   
   for i=1:length(N) 
        mfprintf(fd,"%1.15f\t",N(i));
   end
        mclose(fd);
endfunction
