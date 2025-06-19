## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
## Resubmission

This is a resubmission. In this version, I have addressed the comments from the CRAN team.

## Reviewer Comments

### Regarding DESCRIPTION file formatting

> Please always write package names, software names and API ... in single quotes...

**My Response:**
I have gone through the `DESCRIPTION` file and added single quotes to all package and software names in the `Title` and `Description` fields as requested (e.g., 'Seurat', 'ComplexHeatmap').

> If there are references describing the methods in your package, please add these in the description field...

**My Response:**
I have added the relevant reference describing the methods to the `DESCRIPTION` file in the suggested format: `Gu (2022) <doi:10.1002/imt2.43> and Hao (2024) <doi:10.1038/s41587-023-01767-y>`.

### Regarding invalid URLs

> Found the following (possibly) invalid URLs... Status: 404 Not Found

**My Response:**
I have corrected the broken URLs in the `DESCRIPTION` file (`https://github.com/FanXuRong/SingleCellComplexHeatMap/issues` and `https://github.com/FanXuRong/SingleCellComplexHeatMap/issues` fields) and in the package vignette. The links now correctly point to the active GitHub repository.

### Regarding examples (`\dontrun{}`)

> \dontrun{} should only be used if the example really cannot be executed... Please replace \dontrun with \donttest...

**My Response:**
I have reviewed all examples. I have replaced `\dontrun{}` with `\donttest{}` for examples that take longer than 5 seconds to run, and unwrapped the examples that run quickly. All examples are now available for automated testing.

### Regarding NOTE for new maintainer

Regarding the `NOTE` about a new maintainer, this is my first submission to CRAN.

---

Thank you for your time and feedback.

Best regards,

XueCheng Fang