function patient = TCoMX_MSA(tomoplan, METRICS_LIST, patient)
% TCoMX_MSA Calculates MSA
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

output  = TCoMX_FIND_METRIC(METRICS_LIST, 'MSA');

if ~isempty(output)
    
    if isempty(output.parameters)

	metricname = 'MSA';
    
    Pj = (1:tomoplan.Nleaves)-tomoplan.Nleaves/2+0.5;
    LPS = sum(sinogram,1)*ProjectionTime/TreatmentTime;
    
    % Project the sinogram onto the horizontal axis
    patient.(output.category).(output.subcategory).(metricname) = (sum(LPS.*Pj))/(sum(LPS));

	end
	
end


clear sinogram ProjectionTime TreatmentTime

end
