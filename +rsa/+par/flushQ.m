% flushQ()
%
% This function will delete of the current user's jobs from the CBU queue.
%
% FJ 10/2014

function flushQ()
    myCluster = parcluster('CBU_Cluster');
    delete(myCluster.Jobs);
end
