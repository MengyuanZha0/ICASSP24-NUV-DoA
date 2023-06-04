# DOA Estimation Code

This repository contains the code for DOA (Direction of Arrival) estimation using the NUV-EM (Non-Uniform Variational Expectation-Maximization) algorithm. The code provides implementations for data generation, the NUV-EM algorithm, a baseline algorithm, and the DOA estimation process.

## Table of Contents

- [Introduction](#introduction)
- [Data Generation](#data-generation)
- [NUV-EM Algorithm](#nuv-em-algorithm)
- [Baseline Algorithm](#baseline-algorithm)
- [Estimation](#estimation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Introduction

The DOA estimation problem involves determining the directions of arrival of multiple sources based on measurements from an array of sensors. The NUV-EM algorithm is a novel approach that leverages sparse signal recovery techniques to accurately estimate the DOA in challenging scenarios with limited snapshots and correlated or coherent source signals.

This repository provides the code necessary for data generation, implementing the NUV-EM algorithm, baseline algorithms (MUSIC and Root-MUSIC etc.), and the DOA estimation process.

## Data Generation

The data generation code generates synthetic data with predefined parameters to simulate array measurements. It allows you to customize the number of sources, signal-to-noise ratio (SNR), array number and geometry, and other relevant parameters.

## NUV-EM Algorithm

The NUV-EM algorithm implementation can be found in the `nuv_em.py` file. It includes functions for initializing the algorithm, performing convergence checking of iteration, and construction of final spectrum. The NUV-EM algorithm iteratively estimates the DOA by maximizing the likelihood of the observed measurements.

## Baseline Algorithms

The baseline algorithm code is provided in the `baseline.py` file. It serves as a reference for performance comparison with the NUV-EM algorithm. The baseline algorithm utilizes conventional approaches such as MUSIC or Root-MUSIC.

## Estimation

The estimation code combines the data generation, NUV-EM algorithm, and baseline algorithm to estimate the DOA from the generated measurements. It processes the data, applies the algorithms, and compares the performance to evaluate the effectiveness of the NUV-EM algorithm.

## Usage

To use the code in this repository, follow these steps:

1. Clone the repository: `git clone (https://github.com/mengyuanzhao0812/NUV-EM)`
2. Install the required dependencies: `pip install -r requirements.txt`
3. Customize the data generation parameters in `data_generation.py`
4. Run the estimation code, providing the necessary inputs and configuration options.

Please refer to the specific files and functions for detailed usage instructions and examples.

## Contributing

Contributions to this project are welcome! If you find any issues or have suggestions for improvements, please feel free to open an issue or submit a pull request.

For major contributions, it is recommended to first discuss the proposed changes with the project maintainers through the issue tracker.

## License

This project is licensed under the [MIT License](LICENSE.md).

