function uploadAndCompile
nowStr = datestr(now,30);

% find out if code folder exists
[status,~]  = ...
    system(['ssh dchaimow@guillimin-p2.hpc.mcgill.ca ' ...
    '"[ -d ' fullfile('~','matlab','odcpsfMCMC')  ' ]"']);

% if so, rename the folder
if status == 0
    [status,result]  = ...
        system(['ssh dchaimow@guillimin-p2.hpc.mcgill.ca ' ...
        '"mv ' fullfile('~','matlab','odcpsfMCMC') ...
        ' ' fullfile('~','matlab',['odcpsfMCMC_' nowStr]) '"']);
    if status
        disp(result);
        error('error during code folder renaming');
    end
end

% transfer code folder
[status, result] = system(['scp -r odcpsfMCMC' ...
    ' dchaimow@guillimin-p2.hpc.mcgill.ca:' ...
    fullfile('~','matlab','odcpsfMCMC')]);
if status
    disp(result);
    error('error during code folder transfer');
end

% find out if job folder exists
[status,~]  = ...
    system(['ssh dchaimow@guillimin-p2.hpc.mcgill.ca ' ...
    '"[ -d ' fullfile('~','projectspace','odcpsf','jobs')  ' ]"']);

% if so, rename the folder
if status == 0
    [status,result]  = ...
        system(['ssh dchaimow@guillimin-p2.hpc.mcgill.ca ' ...
        '"mv ' fullfile('~','projectspace','odcpsf','jobs') ...
        ' ' fullfile('~','projectspace','odcpsf',['jobs_' nowStr]) '"']);
    if status
        disp(result);
        error('error during job folder renaming');
    end
end

% compile code
[status,result] = system(['ssh dchaimow@guillimin-p2.hpc.mcgill.ca ' ...
    '"chmod +x ' fullfile('~','matlab','odcpsfMCMC','make.sh') '"']);
if status
    disp(result);
    error('error during compilation');
end

[status,] = system(['ssh dchaimow@guillimin-p2.hpc.mcgill.ca ' ...
    '"' ...
    'source ~/.bash_profile ; ' ...
    'cd ' fullfile('~','matlab','odcpsfMCMC') ' ; '...
    './make.sh' ...
    '"'],'-echo');
if status
    disp(result);
    error('error during compilation');
end
end
