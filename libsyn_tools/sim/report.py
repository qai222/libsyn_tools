from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass
class RunReport:
    summary: dict
    event_log: pd.DataFrame
    shacl_violations: pd.DataFrame
    instance_history: pd.DataFrame

    def write_dir(self, path: str | Path, include_ttl: bool = False) -> None:
        out_dir = Path(path)
        out_dir.mkdir(parents=True, exist_ok=True)

        summary_path = out_dir / "summary.json"
        with summary_path.open("w", encoding="utf-8") as handle:
            json.dump(self.summary, handle, indent=2, sort_keys=True, default=str)

        summary_md_path = out_dir / "summary.md"
        summary_md_path.write_text(self.to_markdown(), encoding="utf-8")

        self.event_log.to_csv(out_dir / "event_log.csv", index=False)

        if include_ttl:
            shacl_df = self.shacl_violations
        else:
            shacl_df = self.shacl_violations.drop(columns=["report_graph_ttl"], errors="ignore")
        shacl_df.to_csv(out_dir / "shacl_violations.csv", index=False)

        self.instance_history.to_csv(out_dir / "instance_history.csv", index=False)

    def to_markdown(self) -> str:
        lines = ["# Run Summary", ""]

        lines.append("## Summary Metrics")
        lines.append("")
        lines.append(f"- makespan: {self.summary.get('makespan')}")
        lines.append("")

        lines.append("## Operation Counts")
        lines.append("")
        op_counts = self.summary.get("operation_counts", {})
        for key in ("start", "end", "abort"):
            lines.append(f"- {key}: {op_counts.get(key, 0)}")
        lines.append("")

        lines.append("## Violations by Shape/Disposition")
        lines.append("")
        lines.append("| shape_iri | disposition | count |")
        lines.append("| --- | --- | --- |")
        for row in self.summary.get("violations_by_shape_disposition", []):
            lines.append(f"| {row.get('shape_iri')} | {row.get('disposition')} | {row.get('count')} |")
        lines.append("")

        lines.append("## Resource Utilization")
        lines.append("")
        utilization = self.summary.get("resource_utilization", {})
        lines.append("### By Pool Type")
        lines.append("")
        for pool_type, count in utilization.get("by_pool_type", {}).items():
            lines.append(f"- {pool_type}: {count}")
        lines.append("")
        lines.append("### By Module")
        lines.append("")
        for module_iri, count in utilization.get("by_module", {}).items():
            lines.append(f"- {module_iri}: {count}")
        lines.append("")

        lines.append("## Remediation Operations")
        lines.append("")
        lines.append(f"- remediation_ops_spawned: {self.summary.get('remediation_ops_spawned', 0)}")
        lines.append("")

        return "\n".join(lines)
