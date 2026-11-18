# Analysis code

Put the code that produces the numbers in `results/` and the figures in the report here. Any language is fine (Python, MATLAB, …).

From the repository root, after downloading both HydRA hulls into `data/`:

```sh
python -m pip install -r requirements.txt
python tools/get_data.py KCS_L0p950
python extra_info/run_all.py
python tools/check_submission.py
```

`msi.py` writes the four official CSV files in `results/`. Research notes, tests, diagnostic figures and the Group 10 answer key live in `extra_info/` and are not part of the course file layout.
