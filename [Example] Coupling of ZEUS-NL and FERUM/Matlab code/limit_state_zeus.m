function g = limit_state_zeus(x1,x2)
	for i=1:length(x1)
        P=x1(i);
        E=x2(i);
        
        % Change input file for Zeus
        fid = fopen('simple.dat'); % open file
        fileContents = textscan(fid,'%s','Delimiter','\n');
        fileContents = fileContents{1};
        fclose(fid); % close file

        % Change modulus of elasticity of material
        matStarts = strmatch('#                                            Material Types',fileContents);
        tmp = regexp(fileContents{matStarts+5}, '\s+', 'split');
        tmp{3} = num2str(E); % column 3 is modulus of elasticity
        tmp = strjoin(tmp);
        fileContents{matStarts+5} = tmp;

        % Change load applied
        loadStarts = strmatch('initial.loads',fileContents);
        tmp = regexp(fileContents{loadStarts+2}, '\s+', 'split');
        tmp{4} = num2str(P); % column 4 is load applied
        tmp = strjoin(tmp);
        fileContents{loadStarts+2} = tmp;

        % write to file 
        fid = fopen('run.dat','wt') ;
        fprintf(fid,'%s\n',fileContents{:});
        fclose(fid) ;

        % Perform structural analysis
        system('RunZEUS.bat');
        count = 0;
        while exist('result.num','file')==0
         pause();
         count = count+1;
         if count == 100
           disp('no file was generated in allowed time');
           break;
         end
        end
        % Read result.num
        rfile = fopen('result.num'); % open file
        resultContents = textscan(rfile,'%s','Delimiter','\n');
        resultContents = resultContents{1};
        fclose(rfile); % close file

        % Determine x - displacement of node 3    
        dispStarts = strmatch('DISPLACEMENTS',resultContents);
        node_ID = 3;
        tmp = regexp(resultContents{dispStarts+1+node_ID}, '\s+', 'split');
        delta_x = str2num(tmp{2});

        system('deleteFILE.bat');

        g(i) = 1.5 - delta_x;
	end
	