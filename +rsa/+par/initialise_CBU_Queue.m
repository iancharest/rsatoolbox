% p = initialise_CBU_Queue(userOptions)
%
% Sets up use of the CBU queueing system.
%
% p: A parpool.
%
% FJ 10-2014
% Jana updated 10-2014

function p = initialise_CBU_Queue(userOptions)

    if userOptions.run_in_parallel || userOptions.run_in_parallel_in_cluster
        
        % Close any existing pool.
        try
            currentPool = gcp('nocreate');
            delete(currentPool);
        catch ex
            prints('parpool initialising ...');
        end
        
        % Can either run in cluster
        if userOptions.run_in_parallel_in_cluster
            cp = cbupool;
            cp.NumWorkers = userOptions.nWorkers;
            cp.SubmitArguments = ['-l walltime=',num2str(userOptions.wallTime), ',mem=' ,num2str(userOptions.memReq),'gb'];      
            if isequal(userOptions.nodesReq , '^N^')
                cp.ResourceTemplate = ['-l nodes=',num2str(userOptions.nodesReq)];    
            else
                cp.ResourceTemplate = ['-l nodes=',num2str(userOptions.nodesReq), ':ppn=' ,num2str(userOptions.proPNode)];    
            end
            p = parpool(cp);
        % Or just on local machine
        else
            p = parpool;
        end  
    end
    
end%function
