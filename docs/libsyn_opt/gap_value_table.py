import pandas as pd


dfs = [
    pd.read_csv("gaps--(A).csv"),
    pd.read_csv("gaps--(B).csv"),
    pd.read_csv("gaps--(C).csv"),
    pd.read_csv("gaps--(D).csv"),
]

df = pd.concat(dfs)
fms_translator = {
    "0": "LAB-1",
    "1": "LAB-2",
}
records = []
for r in df.to_dict("records"):
    name = r["name"]
    percentage_gap = r["percentage gap"]
    solver_status = r["gurobi_status"]
    lib, ntarget, nindex, fms, workshift = name.split("-")
    record = {
        "Library name": ".".join([lib, ntarget, nindex]),
        "Module set": fms_translator[fms],
        "Work shift": bool(int(workshift)),
        "Percentage gap": "N/A" if pd.isna(percentage_gap) else percentage_gap,
        "Optimal": True if solver_status=="Optimal" else False,
    }
    records.append(record)
df = pd.DataFrame.from_records(records)
df.to_csv("gaps_table.csv")