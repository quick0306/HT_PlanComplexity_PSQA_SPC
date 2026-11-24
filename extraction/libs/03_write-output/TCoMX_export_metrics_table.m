function metrics_table = TCoMX_export_metrics_table(patients)
% TCOMX_EXPORT_METRICS_TABLE Flatten computed TCoMX metric structs into a table.
%
% Derived from TCoMX v1.0 (Cavinato et al., J Appl Clin Med Phys 2023;24(1):e13781),
% licensed under the GNU General Public License v3.0 (see extraction/LICENSE).
%
% Local modifications (2026-11-24):
%   - Replaced interactive 'keyboard' call in error handler with 'rethrow(ME)'
%     for non-interactive batch execution.

    metrics_table = [];

    for pidx = 1:numel(patients)

        table_row = [];

        categories = fieldnames(patients(pidx).patient);

        for catidx = 1:numel(categories)

            subcategories = fieldnames(patients(pidx).patient.(categories{catidx}));

            for subcatidx = 1:numel(subcategories)

                try
                    table_row = [table_row struct2table(patients(pidx).patient.(categories{catidx}).(subcategories{subcatidx}))];
                catch ME
                    rethrow(ME);
                end

            end

        end

        metrics_table = [metrics_table; table_row];

    end

end
