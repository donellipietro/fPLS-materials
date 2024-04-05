# fdaPDE Methods Test Bench

## Overview

This repository serves as a test bench specifically designed for evaluating methods related to [fdaPDE](https://fdapde.github.io) (Physics-Informed Spatial and Functional Data Analysis). It provides a collection of utilities and scripts to facilitate the testing process, including model evaluation metrics computation, plot generation, and automation of tests with various parameter configurations.

## Features

- **Model Evaluation Metrics:** Utilities are available for computing various model evaluation metrics (`RMSE`, `IRMSE`, `...`), allowing for comprehensive assessment of fdaPDE methods' performance.
- **Plot Generation:** The repository includes tools for generating plots to visualize the results of the tested methods, aiding in the interpretation and analysis of the experimental outcomes.
- **Test Automation:** Scripts are provided for automating the execution of tests with different parameter configurations. The `run_tests.sh` script facilitates the setup and execution of test suites, streamlining the testing process.

## Usage

### Running Tests

To run tests using the provided utilities, follow these steps::

1. Create a new directory in test containing all the tests scripts.
2. Execute the `run_tests.sh` script, passing the test suite name and test name as arguments. For example:

   ```bash
   ./run_tests.sh my_test_suite my_test
   ```

3. The script initializes the test environment, then iterates over the options files in the specified directory, executing the main test script (main.R) for each file found.
4. After completing the tests, the script performs post-processing tasks, including runtime complexity analysis.

## Repository structure

```bash
.
├── LICENSE
├── Makefile
├── README.md
├── run_tests.sh
├── data
│   └── mesh
│       ├── ...
├── src
│   ├── installation
│   │   ├── install_fdaPDE2.R
│   │   └── install_femR.R
│   └── utils
│       ├── cat.R
│       ├── directories.R
│       ├── meshes.R
│       ├── domain_and_locations.R
│       ├── errors.R
│       ├── results_management.R
│       ├── plots.R
│       └── wrappers.R
├── analysis
└── tests
```

**Files**:

- **LICENSE**: GPL v3 License file specifying the terms and conditions for using the repository.
- **Makefile**: Makefile for automating build tasks or running commands.
- **README.md**: This documentation file providing an overview of the repository and its usage instructions.
- **run_tests.sh**: Script for automating the execution of tests with different parameter configurations.

**Directories**:

- **data**: Directory storing general data files, including mesh data used in tests.
- **src**: Source code directory containing installation scripts and utility functions.
- **tests**: Directory for storing test scripts and related resources.
- **analysis**: Directory containing analysis-related scripts or resources.

## Authorship

This test bench repository is maintained by Pietro Donelli.

## License

This repository is licensed under the GPL v3 License. See the [LICENSE](./LICENSE) file for details.
