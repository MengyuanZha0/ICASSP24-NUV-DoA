# DOA Estimation Code

This repository contains the code for DOA (Direction of Arrival) estimation using the NUV-DoA (Non-Uniform Variational Expectation-Maximization) algorithm. The code provides implementations for data generation, the NUV-EM algorithm, window-sweeping implement and estimation process code, all in Python.

## Table of Contents

- [Introduction](#introduction)
- [Data Generation](#data-generation)
- [NUV-SSR Algorithm](#NUV-SSR)
- [NUV-DoA Algorithm](#NUV-DoA)
- [MUSIC & Root-MUSIC Algorithm](#music)
- [L1-SVD Algorithm](#L1-SVD.mat)
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

## NUV-DoA

high-resolution implementation of NUV-SSR algorithm - NUV-DoA algorithm.

## L1-SVD Algorithm

Implementation of L1-SVD algorithm in MatLab.

## MUSIC & Root-MUSIC Algorithm

Implementation of MUSIC & Root-MUSIC Algorithm.

## Estimation and Interference Canceling

The estimation code combines the data generation, NUV algorithm to estimate the DOA from the generated measurements. We also provides code for interference canceling when faced with multi-sources.


