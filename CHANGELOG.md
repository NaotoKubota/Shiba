# Change log

All notable changes to this Shiba project will be documented in this file.

## [v0.4.0] - 2024-??-??

### Added

### Changed

- Add event type (e.g. `SE`, `FIVE`, `THREE`, `MXE`, `RI`, `MSE`, `AFE`, and `ALE`) to `pos_id` in events and results files.

### Fixed

- Fix a format of `pos_id` of `AFE` and `ALE` events to be consistent between different runs.

## [v0.3.1] - 2024-07-03

### Added

- Add GitHub actions for automatic tag creation and release, and Docker image build and push.

### Changed

### Fixed

## [v0.3] - 2024-07-02

### Added

- New features for identifying multiple skipped exons (MSE), alternative first exon (AFE), and alternative last exon (ALE) ([#4](https://github.com/NaotoKubota/Shiba/issues/4))
  - Up to 5,000 tandemly skipped exons can be identified for MSE, which is reasonable for most datasets and organisms

### Changed

- Unannotated splicing events are now defined when at least one of the associated ***introns*** are not documented in the reference annotation file. Previously, unannotated splicing events were defined when at least one of the associated ***exons*** were not documented in the reference annotation file.
- Update [README.md](https://github.com/NaotoKubota/Shiba/blob/main/README.md) with new badges
- Update [MANUAL.md](https://github.com/NaotoKubota/Shiba/blob/main/doc/MANUAL.md) with new instructions for MSE, AFE, and ALE analyses

### Fixed

- Fixed minor bugs in `SnakeScShiba`

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
