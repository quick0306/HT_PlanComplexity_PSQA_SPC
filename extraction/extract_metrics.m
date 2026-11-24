function extract_metrics(manifest_path, output_path)
% EXTRACT_METRICS Extract plan complexity metrics from TomoTherapy RTPLAN files.
%
% Uses the local TCoMX v1.0 adaptation (GNU GPLv3; Cavinato et al., J Appl Clin
% Med Phys 2023;24(1):e13781) to calculate 10 complexity metrics
% (TTDF, MF, mFLOT, CFNS75, L0NS, CLS, fDISC, centroid, PSTV, MI) and TreatmentTime.
%
% Inputs:
%   manifest_path - CSV file containing columns 'plan_key' and 'rtplan_path'.
%                   Supply 'clinical_dose_per_fraction_gy' for study TTDF.
%                   For a multi-prescription SIB plan, provide the dose per
%                   fraction for the highest prescription level.
%   output_path   - Destination CSV path for extracted complexity metrics table.

root = fileparts(mfilename('fullpath'));
original_directory = pwd;
manifest_path = full_path(manifest_path, original_directory);
output_path = full_path(output_path, original_directory);
manifest_directory = fileparts(manifest_path);

% Read input manifest
input = readtable(manifest_path, FileType="text", Delimiter=",", Encoding="UTF-8", ...
    ReadVariableNames=true, TextType="string", VariableNamingRule="preserve");
keys = string(input.plan_key);
paths = string(input.rtplan_path);

% Ensure clean path restoration upon function exit
original_path = path;
restore_path = onCleanup(@() path(original_path));
restore_directory = onCleanup(@() cd(original_directory));
addpath(genpath(fullfile(root, 'libs')));
cd(root);

% Configure TCoMX metric parameters
archive = TCoMX_READ_METRICS();
config.LOT_statistics.Relative_LOT.mFLOT = archive.LOT_statistics.Relative_LOT.mFLOT;
config.LOT_statistics.Relative_LOT.CFNS.parameters = {'75'};
config.Sinogram.Modulation.PSTV = archive.Sinogram.Modulation.PSTV;
config.Sinogram.Modulation.MI = archive.Sinogram.Modulation.MI;
config.Sinogram.Geometry.LNS.parameters = {'0'};
config.Sinogram.Geometry.CLS.parameters = {'normal'};
config.Sinogram.Geometry.fDISC = archive.Sinogram.Geometry.fDISC;
config.Sinogram.Geometry.centroid = archive.Sinogram.Geometry.centroid;

names = ["TTDF", "MF", "mFLOT", "CFNS75", "L0NS", "CLS", "fDISC", ...
    "centroid", "PSTV", "MI", "TreatmentTime", "TTDF_clinical_s_per_cGy"];
values = NaN(height(input), numel(names));

doses = NaN(height(input), 1);
if ismember('clinical_dose_per_fraction_gy', input.Properties.VariableNames)
    doses = str2double(string(input.clinical_dose_per_fraction_gy));
end

% Process each RTPLAN file
for row = 1:height(input)
    filename = full_path(paths(row), manifest_directory);
    plan = TCoMX_read_RTplan(filename);
    if isempty(plan)
        error('extract_metrics:ReadFailed', ...
            'Failed to parse helical tomotherapy RTPLAN for plan_key "%s" (%s).', ...
            keys(row), filename);
    end
    record.patient = TCoMX_COMPUTE_METRICS(plan, config);
    metrics = TCoMX_export_metrics_table(record);

    % Modulation factor: max leaf open time / mean non-zero leaf open time
    mf = max(plan.sinogram(:)) / mean(nonzeros(plan.sinogram(:)));
    values(row, 2:11) = [mf, metrics{1, names(3:10)}, plan.TreatmentTime];

    % Leave study TTDF missing when the clinical normalization dose is unknown.
    dose = doses(row);
    if isfinite(dose) && dose > 0
        values(row, [1, 12]) = plan.TreatmentTime / (100 * dose);
    end

    if mod(row, 50) == 0 || row == height(input)
        fprintf('Extracted metrics: %d / %d plans.\n', row, height(input));
    end
end

% Write output table
result = [table(keys, VariableNames={'plan_key'}), ...
    array2table(values, VariableNames=cellstr(names))];
writetable(result, output_path);
fprintf('Metrics successfully saved to %s\n', output_path);
end

function filename = full_path(filename, base)
% Convert relative path to absolute path based on base directory
filename = char(filename);
if isempty(regexp(filename, '^(?:[A-Za-z]:[\\/]|[\\/])', 'once'))
    filename = fullfile(base, filename);
end
end
