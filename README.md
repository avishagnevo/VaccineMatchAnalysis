# VaccineMatchAnalysis
# Vaccine Effectivness Analysis on Propencity Scores Balanced Cohort 
> Date: November, 2023. 
## Authors

* Avishag Nevo: avishag.nevo@campus.technion.ac.il


## Project Overview
Welcome to the Epidemiological Impact Analysis repository, a comprehensive reproduction of the published paper [**"BNT162b2 mRNA Covid-19 Vaccine in a Nationwide Mass Vaccination Setting"**](https://www.nejm.org/doi/full/10.1056/nejmoa2101765#article_supplementary_material) by Noa Dagan, MD, et al. led by Professor Yair Goldberg. This statistical project investigates the multifaceted impact of vaccination on infection rates across diverse demographic factors. Leveraging synthetic data generation, sophisticated matching on propensity score algorithms to avoid selection bias, and advanced statistical analysis, this project aims to shed light on the effectiveness of vaccinations in preventing infections over time.

## Project Components

### `main.py`

This script serves as the project's entry point, orchestrating the entire experiment. The `Experiment` class is instantiated, synthetic data is generated, propensity score matching is performed, and statistical analysis is conducted. Execute this script to run the entire experiment seamlessly. You can control the experimant setting by changing the default parametes sent to experiment() function.

### `Experiment.py`

The `Experiment` class is the engine behind data generation for the project. It features methods for data synthesis, propensity score assignment through logistic regression, and compelling visualizations that encapsulate the intricacies of the generated data.


### `Patient.py`

Representing an individual in the experiment, the `Patient` class encapsulates key attributes, including age, vaccination status, infection status, and additional details like vaccination progress and matched counterparts.

### `Matching.py`

The `Matching` class has undergone enhancements, now employing a narrower caliper for refined matching. Visualizations have been elevated to include a scatter plot of propensity scores for matched groups, along with box plots illustrating propensity score intervals. The sequential and dynamic nature of updating treatment and control sets has been optimized for increased efficiency.

### `Analysis.py`

The `Analysis` class undertakes comprehensive statistical analysis, featuring methods for creating detailed tables, conducting Kaplan-Meier survival analysis, Cox proportional hazard modeling, and visualization of standardized differences in means before and after matching. For each person, follow-up ended at the earliest of the following events: occurrence of an outcome event, vaccination (for unvaccinated controls), vaccination of the matched control (for vaccinated persons), or the end of the study period.

## Usage Guide

1. Clone the repository:

   ```bash
   git clone https://github.com/avishagnevo/VaccineMatchAnalysis.git
   cd VaccineMatchAnalysis
   ```

2. Install dependencies for Anaconda/Miniconda:
   
   Project enviorment dependencies listed in `src\enviroment.yml`.
   
   From project root folder, run:
   ```bash
   conda env create -f src/environment.yml -n vma_env
   conda activate vma_env
   ```

3. Run the experiment:

   ```bash
   python src/main.py
   ```

## Dependencies

- `matplotlib`
- `numpy`
- `scikit-learn`
- `lifelines`
- `pandas`
- `tabulate`
- `scipy`


## Visualization Highlights

Immerse yourself in the visual narrative of the project:

- **Histograms**: Uncover patterns in patients demographics.
- **Scatter Plots**: Explore relationships between age, past risk, region risk, infection probabilities and applied propencity scores.
- **Pie Chart**: See the final distribution of infections on vaccinated and control groups.
- **Box Plot**: Delve into matched patients propencity scores proximity distribution.
- **Lollipop Plot**: Check the balancing for each attribute before and after appling matching method.
- **Survival Curves**: Observe Kaplan-Meier survival curves for vaccinated and control groups.
- **Hazard Ratio Dynamics**: Investigate how vaccine effectiveness evolves over time.

## Matching Methodology

The heart of this project lies in the sophisticated matching method employed in `Matching.py`. Key enhancements include:

- **Matching Process**: Every calendar day, all newly vaccinated persons were matched in a 1:1 ratio to unvaccinated controls.  

- **Narrow Caliper**: A refined matching caliper ensures closer matches between treatment and control groups, enhancing the precision of the analysis.

- **Sequential and Dynamic Updates**: Matching was performed using a “rolling cohort” design, the sequential and dynamic updating of treatment and control sets has been optimized for increased efficiency.

- **Visualization**: Visualizations now include scatter plots of propensity scores for matched groups, box plots to understand distanced between matched patients proximity distribution and lollypop chart for the balance achived offering a clearer understanding of the matching process.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
