function [] = set_parpool_cpu_affinity(cpus_per_worker)
% Sets the number of cpu cores that each worker in the parallel pool can use, cpus_per_worker.  If not already started, starts the parallel pool with n_cpu-1 workers. Only works on Windows.  
% If cpus_per_worker*n_works>n_cpus, then some cpus will be shared by multiple workers.
% This function must be re-run whenever the parallel pool is started (i.e. it is not permanent if the parallelpool is shutdown).
% To reset the affinity assignment to default, you can shutdown and restart the parallel pool or simply call set_parpool_cpu_affinity(Inf)
% Only works on Windows because it uses .NET functions

% This function isolates the workers.  Usually performance is about the same, but it can be useful in the following example use cases: 
% Example use cases for this function:
%  (1) Prevent the computer from bogging down by reserving a cpu core: Set the number of workers to 1 less than the number of cpu cores and set cpus_per_worker to 1.  This reserves a core that no workers will use.  This helps when working remotely.
%  (2) If running jobs asynchronously with parfeval() and poll() send(), any one job won't take over the entire computer. 
%  (3) You like the idea of workers having their own core.  Typically this isn't actually faster (performance is usually about the same) and can also be slower. Functions like fft() could otherwise utilize multiple cores, and will be slower when restricted to a single core.

% This code is based on answers by Andrew Janke and others on Stack Exchange. See: https://stackoverflow.com/questions/9778287/set-processor-affinity-for-matlab-engine-windows-7
% which was based on https://stackoverflow.com/questions/2510593/how-can-i-set-processor-affinity-to-a-thread-or-a-task-in-net

    if cpus_per_worker<1
        error('cpus_per_worker<1, cpus_per_worker must be 1 or more.')
    end

    if ispc() % check if windows, othwerwise warn that it doesn't work

        par_pool_started = not(isempty(gcp('nocreate'))); % check if parpool is started

        if par_pool_started

        else % starts a parallel pool with n_cpus-1 workers

            disp('Parallel pool not started yet. Starting parallel pool')
            parpool(feature('numcores')-1); 

        end

        p = gcp(); % get handle for the parallel pool
        poolsize = p.NumWorkers; % get number of workers
        
        parfor i_worker = 1:poolsize % loop through each worker, and set cpu affinity of each
            set_cpu_core(i_worker,cpus_per_worker) % for each worker, set the cpu affinity
        end

    else
        warning('Could not set worker affinity.  This function only works on Windows.')
    end
end
    
function set_cpu_core(i_core,cpus_per_worker)
% This subfunction sets the cpu affinity of each worker in the parallel pool when run using parfor

% This function is edited from Andrew Janke's Stack Exchange Answer: https://stackoverflow.com/questions/9778287/set-processor-affinity-for-matlab-engine-windows-7
% which was based on https://stackoverflow.com/questions/2510593/how-can-i-set-processor-affinity-to-a-thread-or-a-task-in-net

    n_cpus = str2num(getenv('NUMBER_OF_PROCESSORS')); % total number of processors on this computer

    if cpus_per_worker>n_cpus
        cpus_per_worker = n_cpus;
    end
    
    cpu_vec_bin = zeros(1,n_cpus);
    cpu_vec_bin(mod((i_core-1)*cpus_per_worker+(1:cpus_per_worker)-1,n_cpus)+1) = ones(1,cpus_per_worker, 'logical'); % a binary array representing which cpu cores this worker will be assigned
    cpu_dec = sum(2.^(find(cpu_vec_bin)-1)); % convert the binary vector to a decimal
    
    % disp(['Worder #' num2str(i_core) ': ' num2str(cpu_vec_bin) ', ' num2str(cpu_dec)]) % tell the user which cores will be assigned to this worker.
    
    proc = System.Diagnostics.Process.GetCurrentProcess(); % get the current processes
    proc.ProcessorAffinity = System.IntPtr(int32(cpu_dec)); % set affinity mask
    aff = proc.ProcessorAffinity.ToInt32;  % get the affinity mask

    disp(['Worker #' num2str(i_core,'%02.f') ', affinity: ' dec2bin(System.Diagnostics.Process.GetCurrentProcess.ProcessorAffinity.ToInt32, n_cpus)]);

end


