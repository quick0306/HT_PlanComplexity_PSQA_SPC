function patient = TCoMX_SI(tomoplan, METRICS_LIST, patient)
% TCoMX_SI Calculates mSI, sdSI, mdSI.
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
%% Peliminary operations

% Read the sinogram from tomoplan
sinogram = tomoplan.sinogram;

% Read the TreatmentTime from tomoplan
TreatmentTime = tomoplan.TreatmentTime;

% Read the ProjectionTime from tomoplan
ProjectionTime = tomoplan.ProjectionTime;


output  = TCoMX_FIND_METRIC(METRICS_LIST, 'mSI');

if ~isempty(output)
    
    if isempty(output.parameters)

	metricname = 'mSI';
	% Project the sinogram onto the horizontal axis
        patient.(output.category).(output.subcategory).(metricname) = mean(nonzeros(sum(sinogram,1)*ProjectionTime/TreatmentTime)); %*100));

	end
	
end


output  = TCoMX_FIND_METRIC(METRICS_LIST, 'mdSI');

if ~isempty(output)
    
    if isempty(output.parameters)

	metricname = 'mdSI';
	% Compute the SI median
	patient.(output.category).(output.subcategory).(metricname) = median(nonzeros(sum(sinogram,1)*ProjectionTime/TreatmentTime));%*100));

    end
    
end

output  = TCoMX_FIND_METRIC(METRICS_LIST, 'sdSI');

if ~isempty(output)
    
    if isempty(output.parameters)

	metricname = 'sdSI';
	% Compute the SI std 
	patient.(output.category).(output.subcategory).(metricname) = std(nonzeros(sum(sinogram,1)*ProjectionTime/TreatmentTime));%*100));
   
    end
    
end

clear sinogram ProjectionTime TreatmentTime

end
