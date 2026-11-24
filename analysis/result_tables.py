"""Write analysis results as Excel workbooks with readable column widths."""

from pathlib import Path

import pandas as pd
from openpyxl.utils import get_column_letter


def save_table(
    data: dict | list | pd.DataFrame,
    output_dir: Path,
    filename: str,
    *,
    index_label: str | None = None,
) -> None:
    """Save one result table with numeric cells and a frozen header row."""
    table = pd.DataFrame(data)
    if index_label is not None:
        table = table.rename_axis(index_label).reset_index()

    with pd.ExcelWriter(output_dir / f"{filename}.xlsx", engine="openpyxl") as writer:
        table.to_excel(writer, sheet_name="Results", index=False, freeze_panes=(1, 0))
        sheet = writer.sheets["Results"]
        sheet.auto_filter.ref = sheet.dimensions
        for i, name in enumerate(table.columns, start=1):
            values = table[name]
            content_width = int(values.dropna().astype(str).str.len().max()) if values.notna().any() else 0
            width = min(48, max(len(str(name)), content_width, 10) + 2)
            sheet.column_dimensions[get_column_letter(i)].width = width
