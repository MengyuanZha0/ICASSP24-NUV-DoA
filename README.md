# DOA Estimation Code

This repository contains the code for DOA (Direction of Arrival) estimation using the NUV-DoA (Non-Uniform Variational Expectation-Maximization) algorithm. The code provides implementations for data generation, the NUV-EM algorithm, window-sweeping implement and estimation process code, all in Python.

## Table of Contents

- [Introduction](#introduction)
- [Data Generation](#data-generation)
- [NUV-SSR Algorithm](#NUV-SSR)
- [NUV-DoA Algorithm](#NUV-DoA)
- [Estimation](#estimation)

## Introduction

The DOA estimation problem involves determining the directions of arrival of multiple sources based on measurements from an array of sensors. The NUV-EM algorithm is a novel approach that leverages sparse signal recovery techniques to accurately estimate the DOA in challenging scenarios with limited snapshots and correlated or coherent source signals.

This repository provides the code necessary for data generation, implementing the NUV-EM algorithm, baseline algorithms (MUSIC and Root-MUSIC etc.), and the DOA estimation process.

## Data Generation

Generating data for this experiment consists of the following steps:
1. Generate the true DoAs, and generate its corresponding steering matrix.
2. Generate noncoherent source signals.
3. Define noise level, and generate the observations

## NUV-SSR

This file presents the original NUV-SSR algorithm benchmark.

## Window-sweeping Implementation

high-resolution implementation of NUV-SSR algorithm - NUV-DoA algorithm.

## Estimation and Interference Canceling

The estimation code combines the data generation, NUV algorithm to estimate the DOA from the generated measurements. We also provides code for interference canceling when faced with multi-sources.


