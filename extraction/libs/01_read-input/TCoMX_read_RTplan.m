function tomoplan = TCoMX_read_RTplan(filename, varargin)
% TCOMX_READ_RTPLAN Read a helical tomotherapy DICOM RTPLAN file.
%
% Derived from TCoMX v1.0 (Cavinato et al., J Appl Clin Med Phys 2023;24(1):e13781),
% licensed under the GNU General Public License v3.0 (see extraction/LICENSE).
%
% Local modifications (2026-11-24):
%   - Carry forward dynamic jaw width and sinogram leaf-open times when DICOM
%     control points omit repeated values.
%   - Tolerate missing optional or de-identified patient/operator metadata tags.

try
    % Read DICOM RT
    info = dicominfo(filename);
    if isfield(info, 'Modality') && strcmpi(info.Modality, 'RTPLAN')
        % Define the plan name, Add HA prefix for Hi-Art plans
        plan_name = get_field_or_default(info, 'RTPlanName', '');
        tomoplan.TPS = 'TOMO Legacy';
        tomoplan.name = ['HA_' char(plan_name)];

        %% GENERAL INFORMATION (optional for de-identified DICOM RTPLAN files)

        tomoplan.OriginalPlanName = plan_name;
        tomoplan.PatientID = get_field_or_default(info, 'PatientID', '');
        tomoplan.PatientName = get_field_or_default(info, 'PatientName', struct('FamilyName', ''));
        tomoplan.PatientBirthDate = get_field_or_default(info, 'PatientBirthDate', '');
        tomoplan.PatientSex = get_field_or_default(info, 'PatientSex', '');
        tomoplan.RTPlanDate = get_field_or_default(info, 'RTPlanDate', '');
        reviewer = get_field_or_default(info, 'ReviewerName', struct('FamilyName', ''));
        tomoplan.Physicist = get_field_or_default(reviewer, 'FamilyName', '');
        operator = get_field_or_default(info, 'OperatorsName', struct('FamilyName', ''));
        tomoplan.Physician = get_field_or_default(operator, 'FamilyName', '');

        %% IMPORT ATTRIBUTES THAT ARE EQUAL IN RAYSTATION AND HI-ART

        % Import plan manufacturer name
        tomoplan.ManufacturerModelName = get_field_or_default(info, 'ManufacturerModelName', '');

        % Read the number of CP
        tomoplan.NCP = info.BeamSequence.Item_1.NumberOfControlPoints;

        % Define the rotating direction
        tomoplan.RotationDirection = info.BeamSequence.Item_1.ControlPointSequence.Item_1.GantryRotationDirection;

        % Allocate the vector with the gantry angles
        tomoplan.GantryAngle = zeros(tomoplan.NCP, 1);

        % Import number of leaves
        tomoplan.Nleaves = 64;

        % Import information at each CP
        tomoplan.sinogram = zeros(tomoplan.NCP, tomoplan.Nleaves);

        % Number of Projections per rotation
        tomoplan.NProjectionsPeriod = 51;

        % Number of Rotations
        tomoplan.NRotations = (tomoplan.NCP - 1) / tomoplan.NProjectionsPeriod;

        % Import the Gantry Period [s]
        tomoplan.GantryPeriod = str2double(strsplit(char(info.BeamSequence.Item_1.Private_300d_1040), '\'));

        % Compute the ProjectionTime from the Gantry Period [s]
        tomoplan.ProjectionTime = tomoplan.GantryPeriod / tomoplan.NProjectionsPeriod;

        % Compute the TreatmentTime [s]
        tomoplan.TreatmentTime = tomoplan.NRotations * tomoplan.GantryPeriod;

        % Import the CouchSpeed [mm/s]
        tomoplan.CouchSpeed = str2double(strsplit(char(info.BeamSequence.Item_1.Private_300d_1080), '\'));

        % Compute the CouchTranslation [mm]
        tomoplan.CouchTranslation = tomoplan.CouchSpeed * tomoplan.TreatmentTime;

        % Import the Pitch
        tomoplan.NominalPitch = str2double(strsplit(char(info.BeamSequence.Item_1.Private_300d_1060), '\'));

        % Compute the Field Width [mm]
        tomoplan.NominalFieldWidth = tomoplan.CouchTranslation / tomoplan.NRotations / tomoplan.NominalPitch;

        % Compute the TargetLength [mm]
        tomoplan.TargetLength = tomoplan.CouchTranslation - tomoplan.NominalFieldWidth;

        % Import BeamMeterset (predicted duration of the treatment)
        tomoplan.BeamMeterset = info.FractionGroupSequence.Item_1.ReferencedBeamSequence.Item_1.BeamMeterset;

        %% IMPORT ATTRIBUTES THAT ARE DIFFERENT IN HI-ART
        for cp = 1:tomoplan.NCP
            cp_item = info.BeamSequence.Item_1.ControlPointSequence.(['Item_' num2str(cp)]);

            % Import the GantryAngle at the current CP
            tomoplan.GantryAngle(cp, 1) = cp_item.GantryAngle;

            % Import the CouchPosition at the current CP
            tomoplan.CouchPosition(cp, 1) = cp_item.IsocenterPosition(end);

            if isfield(info.BeamSequence.Item_1, 'BeamType') && strcmpi(info.BeamSequence.Item_1.BeamType, 'DYNAMIC')
                % In dynamic-jaw mode, DICOM RTPLAN omits unchanged jaw positions at subsequent CPs.
                try
                    if cp == 1
                        jaw_pos = cp_item.BeamLimitingDevicePositionSequence.Item_2.LeafJawPositions;
                    else
                        jaw_pos = cp_item.BeamLimitingDevicePositionSequence.Item_1.LeafJawPositions;
                    end
                    tomoplan.inFieldWidth(cp, 1) = jaw_pos(2) - jaw_pos(1);
                catch
                    if cp > 1
                        tomoplan.inFieldWidth(cp, 1) = tomoplan.inFieldWidth(cp - 1, 1);
                    else
                        tomoplan.inFieldWidth(cp, 1) = tomoplan.NominalFieldWidth;
                    end
                end
            end

            try
                % Import the plan sinogram at the current CP (for all leaves).
                % DICOM RTPLAN omits repeated MLC leaf open times when unchanged from the previous CP.
                tmp = char(cp_item.Private_300d_10a7');
                if isempty(tmp)
                    if cp > 1
                        tomoplan.sinogram(cp, :) = tomoplan.sinogram(cp - 1, :);
                    end
                else
                    tmp = strsplit(tmp, '\');
                    tomoplan.sinogram(cp, :) = str2double(tmp);
                end
            catch
                % Leave zero row when private sinogram tag is absent for a CP.
            end
        end

        try
            % Import the fraction dose (works for RayStation)
            tomoplan.DoseFraction = info.FractionGroupSequence.Item_1.ReferencedDoseReferenceSequence.Item_1.TargetPrescriptionDose * 100 / info.FractionGroupSequence.Item_1.NumberOfFractionsPlanned;
            tomoplan.TTDF = tomoplan.TreatmentTime / tomoplan.DoseFraction;
        catch
            tomoplan.TTDF = tomoplan.TreatmentTime / 200;
        end
    else
        tomoplan = [];
    end

catch
    tomoplan = [];
end

end

function value = get_field_or_default(s, field_name, default_value)
if isstruct(s) && isfield(s, field_name)
    value = s.(field_name);
else
    value = default_value;
end
end
