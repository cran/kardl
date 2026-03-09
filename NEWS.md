# kardl 1.1.0

## Changes in this version

- **Standard classes:** All tests and estimation routines have been rewritten to use standard R classes such as `htest`, `anova`, and `lm`. This improves compatibility with generic methods and downstream analyses.

- **Multipliers & bootstrap:** Added functionality to compute dynamic multipliers along with a bootstrap-based inference method, allowing more robust uncertainty quantification.

- **ECM method updated:** The error correction model (ECM) method has been modified to enhance estimation accuracy and better integrate with the new class structure.

- **Internal improvements:** Minor internal code improvements and documentation updates to support new features and maintain CRAN compliance.

- **Windows compatibility:** Fixed an issue where the `summary()` output was not fully compatible with Windows in earlier versions.

- **Documentation:** Updated documentation to reflect changes in function names and new features, ensuring clarity for users.

- **Examples:** Updated examples in the documentation to demonstrate the new features and changes in function usage.

- **Testing:** Added new unit tests to cover the new features and ensure the robustness of the package.

- **Performance:** Improved performance of the estimation routines, particularly for larger datasets, through optimized code and better use of R's vectorized operations.

- **Error handling:** Enhanced error handling to provide more informative messages when users encounter issues with input data or function usage.

- **Vignettes:** Updated vignettes to include examples of the new features and to provide a comprehensive overview of the package's capabilities.
