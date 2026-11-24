function patient = TCoMX_nOC(tomoplan, METRICS_LIST, patient)
% TCoMX_nOC Calculates nOC
%
%   INPUT ARGUMENTS:
%      tomoplan: structure array containing information about a
%                TomoTherapy® plan as returned by the function TCoMX_read_plan
%                of TCoMX.
%       METRICS_LIST: structure array containing all the metrics to be
%                     computed according to the METRICS.in file.
%       patient: structure containing the results of the computation
%
%   OUTPUT ARGUMENTS:
%       patient: structure array containing the metrics computed according
%                to the METRICS.in file. The sub-fields are organized in
%                categories and sub-categories accordingly to the "TomoTherapy® 
%                Complexity Metrics Reference Guide"
%
%% 
output  = TCoMX_FIND_METRIC(METRICS_LIST, 'nOC');

if ~isempty(output)

    if isempty(output.parameters)
        sinogram = tomoplan.sinogram;

        sinogram(sinogram == 1)  = -1;
        sinogram(sinogram > 0)   = 1;
        sinogram(sinogram == -1) = 2;

        leafmov = zeros(1, tomoplan.Nleaves);

        for nl=1:tomoplan.Nleaves
            for cp = 1:tomoplan.NCP-1
                if sinogram(cp, nl) == 0
                    if sinogram(cp+1, nl) == 2
                        leafmov(nl) = leafmov(nl)+1;  
                    end
                end

                if sinogram(cp,nl) == 1
                    leafmov(nl) = leafmov(nl)+2;
                    if sinogram(cp+1, nl) == 2
                        leafmov(nl) = leafmov(nl)+1;
                    end
                end

                if sinogram(cp, nl) == 2
                    if sinogram(cp+1, nl) == 0
                        leafmov(nl) = leafmov(nl)+1;
                    elseif sinogram(cp+1, nl) == 1
                        leafmov(nl) = leafmov(nl)+1;
                    end
                end    
            end
        end

        patient.(output.category).(output.subcategory).(output.metric)  = mean(leafmov)/(tomoplan.NCP-1);

    end
end

end