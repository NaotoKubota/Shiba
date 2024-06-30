# Change log

All notable changes to this Shiba project will be documented in this file.

## [v0.2.2] - 2024-06-30

### Added

### Changed

### Fixed

- Fix `numpy` dependency in the Docker image to `1.26.4`

## [v0.2.1] - 2024-06-29

### Added

### Changed

### Fixed

- Fix the bug in `SnakeShiba` ([#7](https://github.com/NaotoKubota/Shiba/issues/7))

## [v0.2] - 2024-06-17

### Added

- `pca.py`: Principal Component Analysis (PCA) module ([#1](https://github.com/NaotoKubota/Shiba/issues/1))
  - Perform PCA for top 3,000 highly variable genes or splicing events
    - When PSI matrix has less than 6,000 rows without NaN values, drop rows with more than 50% NaN values followed by KNN imputation (n_neighbors=5), and then get top 3,000 highly variable splicing events
- This `CHANGELOG.md` file to keep track of changes

### Changed

- Generate PSI matrix (`PSI_matrix_sample.txt` and `PSI_matrix_group.txt`) when differential splicing analysis is performed
- Update the style of `summary.html`
- Update [MANUAL.md](https://github.com/NaotoKubota/Shiba/blob/main/doc/MANUAL.md) with new instructions
- Update [README.md](https://github.com/NaotoKubota/Shiba/blob/main/README.md) with citation information

### Fixed

## [v0.1] - 2024-05-17

### Added

- Initial release

### Changed

### Fixed
