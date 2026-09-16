# Environment Notes

Known final validation environment:

- R: 4.5.2
- Python: project-specific virtual environment used for Ag3 acquisition
- Windows PowerShell used for the acquisition wrapper

Core R packages used by the final validation scripts:

```text
data.table
ggplot2
readr
```

The final Step-07 script also writes `sessionInfo.txt` to the output directory.

For the public release, consider generating:

```r
sessionInfo()
installed.packages()[, c("Package", "Version")]
```

and recording the Python environment with:

```powershell
python --version
python -m pip freeze
```

Do not commit environment folders or credentials.
