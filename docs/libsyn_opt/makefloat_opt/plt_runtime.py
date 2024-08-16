# def plot_runtime_ax(runs: list[SchedulerRun], subfig_title: str, ax: plt.Axes, remove_legend):
#     df = pd.DataFrame.from_records([r.model_dump() for r in runs])
#     p = sns.stripplot(
#         df, x="n_target", y="runtime_milp", ax=ax, log_scale=True,
#         # dodge = True,
#         # orient = "h",
#         alpha=.5, s=10, jitter=False, linewidth=1, marker="x"
#
#     )
#     ax.axhline(y=3600, label="1 hr", color="red", ls=":")
#     # ax.axhline(y=10800, label="3 hr", color="red")
#     ax.set_ylabel("MILP runtime")
#     ax.set_xlabel("Number of target chemicals")
#     ax.legend()
#     if remove_legend:
#         try:
#             ax.get_legend().remove()
#         except AttributeError:
#             pass
#     # ax.xaxis.set_major_locator(tck.MultipleLocator())
#     # ax.yaxis.set_major_locator(tck.MultipleLocator(base=1))
#     # ax.set_xlim([0, 20])
#     # ax.set_title(subfig_title, fontsize='large', loc='center')
#     ax.set_title(subfig_title, loc='left')
#
#
# def plot_runtime(runs: list[SchedulerRun], figname: str):
#     runs_fda_ws0, runs_fda_ws1, runs_vs_ws0, runs_vs_ws1 = split_runs(runs)
#     # runs_fda = runs_fda_ws0 + runs_fda_ws1
#     # runs_vs = runs_vs_ws0 + runs_vs_ws1
#     fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(8, 8), sharey=True, sharex=True)
#     plot_runtime_ax(runs_fda_ws0, subfig_title="(A) FDA", ax=ax1, remove_legend=False)
#     plot_runtime_ax(runs_fda_ws1, subfig_title="(B) FDA + work shifts", ax=ax2, remove_legend=True)
#     plot_runtime_ax(runs_vs_ws0, subfig_title="(C) VS", ax=ax3, remove_legend=True)
#     plot_runtime_ax(runs_vs_ws1, subfig_title="(D) VS + work shifts", ax=ax4, remove_legend=True)
#     fig.tight_layout()
#     fig.savefig(figname, dpi=600)
