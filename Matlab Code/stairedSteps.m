function u = stairedSteps(samplingTime,ns,amplitude,clock)

%A function for creating ZMP reference (control)
    % samplingTime --> e.g. t = 0:samplingTime:Tf
    % ns           --> number of steps
    % amplitude    --> vector of amplitudes of all steps
    % clock        --> vector of times at which the steps occur
    
    ratio = ns/samplingTime;
    value = 0;
    count = 1;
    u = zeros(1,ratio+1);
    for i=1:ratio+1
        for j=count:count
            if mod(i,clock(j)*ratio/ns) == 0
                value = 0;
                if count < ns
                    count = count + 1;
                end
                for h=1:j
                    value = value + amplitude(h);                     
                end    
            end 
        end
        u(i) = value;
    end
end

