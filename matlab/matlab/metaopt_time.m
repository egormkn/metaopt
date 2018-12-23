function varargout=metaopt_time(name, timeout, varargin)
% default function call with timeout to call the metaopt lib
% timeout has to be given as a string ending with s for seconds, 
% m for minutes, h for hours or d for days

	% this function uses files to transport the parameters to the metaopt lib
	% if you don't like the file names or their location, you can change them here
	% for example, if you have access to a RAM-disk, it may be advantagous to write the files there, 
	% because then you omit writing to harddisk, which is an enormous performance-bottleneck
    % the paths are realtive to the location of this script
	infile_b = 'metaopt_input_'; % file name beginning of the file used for writing input parameters
	outfile_b = 'metaopt_output_'; % file name beginning of the file used for writing output parameters
	
    r = random('unid',1000000);
    infile = [infile_b num2str(r) '.mat'];
    outfile = [outfile_b num2str(r) '.mat'];
	
	% these are the environment variables that are used to tell the metaopt lib 
	% where to look for the input arguments
	% and where to write the output
	% env_file_name_input = 'metaopt_input_file';
	% env_file_name_output = 'metaopt_output_file';
	
	% we will also use an environment variable to tell metaopt which function is actually called
	% env_function = 'metaopt_function';
	
	% now we can start the actual work
    
    metaopt_fullpath = mfilename('fullpath');
    metaopt_name = mfilename;
    metaopt_path=metaopt_fullpath(1:length(metaopt_fullpath)-length(metaopt_name));
	
    try
        % write the input parameters
        save([metaopt_path infile], 'varargin');

        outargs = nargout;
        if outargs == 0
            outargs = 1;
        end

        syscall = ['LD_LIBRARY_PATH="" && ' 'timeout ' timeout ' ' metaopt_path '../bin/metaopt ' metaopt_path infile ' ' metaopt_path outfile ' ' name ' ' num2str(outargs)]

        [status, output] = system(syscall,'-echo');

        if status ~= 0
            disp(output);
            disp('metaopt call failed');
            if status == 124
                disp('reason: timeout');
                ME = MException('metaopt:error', sprintf('aborted computation due to timeout'));
                throw(ME);
            else
            ME = MException('metaopt:error', sprintf('error code: %d',status));
            throw(ME);
            end
        else
            out = load([metaopt_path outfile],'varargout');
            varargout = out.varargout;
        end
        % delete the files used for transporting parameters
        delete([metaopt_path infile]);
        delete([metaopt_path outfile]);
    catch ME
        % make sure the files get deleted (MATLAB doesn't know finally)
        delete([metaopt_path infile]);
        delete([metaopt_path outfile]);
        rethrow(ME);
    end
end