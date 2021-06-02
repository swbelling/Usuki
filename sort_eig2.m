function [sortedvals,sortedvecs,numforwardprop] = sort_eig2(vectors,values)
    Nx = size(values,1);
    %Ny = size(values,2);
    logvalues = log(values);
    backev_vals = double.empty;
    forwardprop_vals = double.empty;
    forwardprop_vecs = double.empty;
    backev_vecs = double.empty;
    forwardev_vals = double.empty;
    forwardev_vecs = double.empty;
    backprop_vals = double.empty;
    backprop_vecs = double.empty;
    for i = 1:Nx
        if abs(abs(values(i,i)) -1) > 1e-6
            if abs(values(i,i)) > 1 
                backev_vals = [backev_vals,values(i,i)];
                backev_vecs = [backev_vecs,vectors(:,i)];
            else
                forwardev_vals = [forwardev_vals,values(i,i)];
                forwardev_vecs = [forwardev_vecs,vectors(:,i)];
            end
        else
            if imag(values(i,i)) < 0 
                backprop_vals = [backprop_vals,values(i,i)];
                backprop_vecs = [backprop_vecs,vectors(:,i)];
            else
                forwardprop_vals = [forwardprop_vals,values(i,i)];
                forwardprop_vecs = [forwardprop_vecs,vectors(:,i)];
            end
        end
    end
    numforwardprop = size(forwardprop_vals,2);
    sortedvals = [forwardprop_vals,forwardev_vals,backprop_vals,backev_vals];
    sortedvecs = [forwardprop_vecs,forwardev_vecs,backprop_vecs,backev_vecs];