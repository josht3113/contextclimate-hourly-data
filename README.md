# contextclimate-hourly-data

Hourly ASOS archive for the [ContextClimate](https://contextclimate.io) Streak Tracker: one compact
file per station per decade for every station in the site's shared roster, so the tool never has to pull
decades of observations from IEM on each visit. The tool fetches only the years after `through` live.

- Served at `https://josht3113.github.io/contextclimate-hourly-data/` (`index.json`, `hourly/<ICAO>/<decade>.bin.gz`).
- Source: Iowa Environmental Mesonet ASOS archive (routine and special reports).
- Builder: `scripts/build_hourly.py` (format and hourly rules are documented at the top of the file).
- Workflow: `Build hourly archive` — run once with mode **full**; it then updates itself every January.
  If a job stops at GitHub's time limit, choose **Re-run failed jobs** and it resumes.
- To add stations after changing the roster, run mode **incremental** (new stations get their full history)
  or mode **stations** with a list.
