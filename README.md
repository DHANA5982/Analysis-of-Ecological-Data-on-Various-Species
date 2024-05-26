# Analysis of Ecological Data on Various Species

## Overview
This repository contains the analysis of ecological data focusing on five species: Bees, Hoverflies, Isopods, Ladybirds, and Grasshoppers & Crickets. The analysis was conducted using various statistical methods and machine learning techniques to explore the relationships and distributions of these species across different land classes in Wales and Scotland over two different periods.

## Table of Contents
- [Introduction](#introduction)
- [Data](#data)
- [Methods](#methods)
- [Results](#results)
  - [Summary Statistics](#summary-statistics)
  - [Correlation Analysis](#correlation-analysis)
  - [Boxplots](#boxplots)
  - [Hypothesis Tests](#hypothesis-tests)
  - [Linear Regression](#linear-regression)
- [Conclusion](#conclusion)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Introduction
This project investigates the ecological data of five species, analyzing their distributions and relationships across various land classes. The study includes summary statistics, winsorised means, correlation analysis, and hypothesis testing, along with linear regression models to predict the richness of these species.

## Data
The data consists of measurements for Bees, Hoverflies, Isopods, Ladybirds, and Grasshoppers & Crickets across different land classes in Wales and specific locations in Scotland. The dataset includes:
- Min, 1st Quartile, Median, Mean, 3rd Quartile, Max, and Winsorised Mean values for each species.
- Correlation coefficients between the species.
- Boxplot comparisons across two periods.
- Results from hypothesis tests (Kolmogorov-Smirnov Test, t-test).
- Linear regression models (simple and multiple).

## Methods
The analysis was conducted using the following methods:
- **Summary Statistics**: Calculated for each species to understand their distributions.
- **Winsorised Mean**: Used to handle outliers by replacing the extreme values.
- **Correlation Analysis**: Assessed the linear relationship between different species.
- **Boxplots**: Visualized the central tendency and variability of species across two periods.
- **Hypothesis Tests**: Conducted Kolmogorov-Smirnov test and t-test to compare distributions and means.
- **Linear Regression**: Applied simple and multiple linear regression to predict species richness based on ecological factors.

## Results

### Summary Statistics
- **Bees**: Min: 0.10, 1st Quartile: 0.37, Median: 0.63, Mean: 0.68, 3rd Quartile: 0.91, Max: 3.31, Winsorised Mean: 0.64
- **Hoverflies**: Min: 0.42, 1st Quartile: 0.66, Median: 0.75, Mean: 0.75, 3rd Quartile: 0.83, Max: 1.00, Winsorised Mean: 0.76
- **Isopods**: Min: 0.18, 1st Quartile: 0.46, Median: 0.69, Mean: 0.65, 3rd Quartile: 0.86, Max: 1.00, Winsorised Mean: 0.65
- **Ladybirds**: Min: 0.12, 1st Quartile: 0.60, Median: 0.75, Mean: 0.74, 3rd Quartile: 0.87, Max: 1.38, Winsorised Mean: 0.74
- **Grasshoppers & Crickets**: Min: 0.07, 1st Quartile: 0.43, Median: 0.56, Mean: 0.59, 3rd Quartile: 0.76, Max: 1.32, Winsorised Mean: 0.59

### Correlation Analysis
- Bees and Hoverflies: 0.06
- Bees and Isopods: -0.04
- Bees and Ladybirds: 0.37
- Bees and Grasshoppers & Crickets: 0.24
- Hoverflies and Isopods: 0.37
- Hoverflies and Ladybirds: 0.41
- Hoverflies and Grasshoppers & Crickets: 0.43
- Isopods and Ladybirds: 0.07
- Isopods and Grasshoppers & Crickets: 0.57
- Ladybirds and Grasshoppers & Crickets: 0.25

### Boxplots
- Compared the distribution of Ladybirds across two periods, showing a right-skew in the earlier period and asymmetric distribution in the later period.

### Hypothesis Tests
- **Kolmogorov-Smirnov Test**: Indicated that the two data sets come from different distributions (p-value < 0.05).
- **t-test**: Suggested that the mean values of BD5 change and BD11 change are equal (p-value > 0.05).

### Linear Regression
- **Simple Linear Regression**: Used Carabids as the response variable and ecological status 5 as the predictor variable. Found significant slope but low R-squared value.
- **Multiple Linear Regression**: Initially included all five species, then refined the model by dropping Ladybirds and adding interaction terms. The final model showed improved R-squared and AIC values.

## Conclusion
The analysis revealed significant insights into the distributions and relationships of the five species. The final multiple linear regression model provided the best fit, indicating potential predictors for species richness.

