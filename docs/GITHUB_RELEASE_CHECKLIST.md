# GitHub Release Checklist

Before making the repository public:

- [ ] Verify target-lock SHA-256
- [ ] Confirm Ag3 did not alter discovery ranks or thresholds
- [ ] Remove credentials, tokens, local usernames, and absolute machine paths
- [ ] Confirm no restricted/raw genomic data are accidentally committed
- [ ] Include README and methods documentation
- [ ] Include `sessionInfo.txt`
- [ ] Record Python package versions
- [ ] Confirm all scripts use relative/project paths where possible
- [ ] Add a software license after author agreement
- [ ] Complete `CITATION.cff`
- [ ] Add repository URL to `CITATION.cff`
- [ ] Tag the first archival release
- [ ] Consider Zenodo linkage for DOI assignment
- [ ] Archive compact processed data needed to reproduce tables/figures
- [ ] Keep large public raw data external and reacquire via scripts
