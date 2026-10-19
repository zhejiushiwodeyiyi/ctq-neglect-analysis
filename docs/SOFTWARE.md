# Software dependencies

The historical audit records in the local research package reported Python 3.12.14, NumPy 2.3.5, pandas 2.2.3, SciPy 1.18.1, statsmodels 0.15.0, Patsy 1.0.3 and scikit-learn 1.9.1. The historical version of threadpoolctl was not recorded. These values reproduce the supplied audit metadata; they were not checked against package-index availability or installed as part of this GitHub preparation.

Current figure preparation also uses matplotlib, Pillow and PyMuPDF. `code/runtime_versions.json` preserves the supplied version record with machine-specific executable and module paths removed. The statistical environment and document/figure preparation environment are separate records.

`requirements.txt` is an unpinned direct-dependency list. Installing it does not establish that the resulting versions will exactly reproduce archived numerical outputs. Record the actual environment and compare outputs before claiming numerical reproduction. No dependency installation or statistical model fitting was performed when this repository was prepared.

The data-free repository check uses only the Python standard library and can run before any dependency installation.
