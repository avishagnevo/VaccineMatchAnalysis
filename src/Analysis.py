import itertools
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
import pandas as pd
import numpy as np
from tabulate import tabulate
from scipy.stats import pearsonr


class Analysis:
    def __init__(self, experiment, matcher):
        self.matched_treatment = matcher.matched_treatment
        self.matched_control = matcher.matched_control
        self.infected_matched_treatment = None
        self.infected_matched_control = None
        self.followup_duration = experiment.followup_duration

        self.survival_matched_treatment = [patient.get_survival_data() for patient in
                                           self.matched_treatment]
        self.survival_matched_control = [patient.get_survival_data() for patient in
                                         self.matched_control]
        self.survival_matched_population = self.survival_matched_treatment + self.survival_matched_control

        self.event_time_treatment = [patient_survival[0] for patient_survival in
                                     self.survival_matched_treatment]
        self.event_time_control = [patient_survival[0] for patient_survival in
                                   self.survival_matched_control]

        self.censored_treatment = [patient_survival[1] for patient_survival in
                                   self.survival_matched_treatment]
        self.censored_control = [patient_survival[1] for patient_survival in
                                 self.survival_matched_control]

    def kaplan_meier(self):
        # Fit Kaplan-Meier survival curves
        kmf_vaccinated = KaplanMeierFitter()
        kmf_vaccinated.fit(self.event_time_treatment,
                           event_observed=self.censored_treatment)

        kmf_unvaccinated = KaplanMeierFitter()
        kmf_unvaccinated.fit(self.event_time_control,
                             event_observed=self.censored_control)

        '''
        kmf_21_v = kmf_vaccinated.cumulative_density_at_times(21)
        kmf_42_v = kmf_vaccinated.cumulative_density_at_times(42)
        kmf_21_u = kmf_unvaccinated.cumulative_density_at_times(21)
        kmf_42_u = kmf_unvaccinated.cumulative_density_at_times(42)
        print(kmf_21_v, kmf_42_v, kmf_21_u, kmf_42_u)
        '''

        # Plot survival curves
        plt.figure(figsize=(12, 6))
        kmf_vaccinated.plot_survival_function(label='Treatment', color='steelblue')
        kmf_unvaccinated.plot_survival_function(label='Control', color='indianred')
        plt.title('Survival Curves for Vaccinated and Unvaccinated Groups')
        plt.xlabel('Days')
        plt.ylabel('Survival Probability')
        plt.legend()
        plt.show()

    def cox_proportional_hazard(self, include_first_period=False):
        covariates = ['age', 'past', 'region', 'vaccinated']
        data = [[censored, time] + [
            getattr(patient, covariate) for covariate in covariates] for
                time, censored, patient in
                self.survival_matched_population]

        matched_data_df = pd.DataFrame(data,
                                       columns=['censored',
                                                'event_time'] + covariates)

        cox_model = CoxPHFitter()  # Calculate hazard ratios using Cox Proportional-Hazards model
        cox_model.fit(matched_data_df, duration_col='event_time',
                      event_col='censored')
        cox_model.print_summary()

        plt.figure(figsize=(12, 6))
        cox_model.plot(hazard_ratios=False, c='black')
        plt.title('Coefficients Log Hazard Ratios with 95%CI')
        plt.xlabel('log(Hazard Ratio)')
        plt.ylabel('Coefficients')
        plt.show()

        maximal_observed_event_time = matched_data_df['event_time'].max()

        followup_times = list(
            range(maximal_observed_event_time)) if include_first_period else list(
            range(14, maximal_observed_event_time))

        hazard_df = cox_model.compute_followup_hazard_ratios(matched_data_df,
                                                             followup_times)
        hazard_df['Vaccine Effectiveness'] = 1 - hazard_df['vaccinated']
        print(hazard_df.round(5))

        # Custom plot for hazard ratios over time
        plt.figure(figsize=(12, 6))
        plt.plot(followup_times, hazard_df['vaccinated'],
                 label='Vaccine Hazard Ratio', marker='o', c='black')
        plt.axhline(1, linestyle='--', color='gray', label='Reference (No Effect)')
        plt.title('Hazard Ratio Over Time (T/C)')
        plt.xlabel('Time')
        plt.ylabel('Hazard Ratio')
        plt.legend()
        plt.show()

        plt.figure(figsize=(12, 6))
        plt.plot(followup_times, hazard_df['Vaccine Effectiveness'], marker='o',
                 c='black')
        plt.title('Vaccine Effectiveness Over Time')
        plt.xlabel('Time')
        plt.ylabel('1 - Vaccine Hazard Ratio')
        plt.show()

    def create_table_1(self):  # sourcery skip: merge-list-appends-into-extend
        covariates = ['age', 'past', 'region']
        characteristics = ['mean', 'std', 'range']
        values_treated_matched = {}
        values_control_matched = {}

        for covariant in covariates:
            values_treated_matched[(covariant, 'values')] = [getattr(patient, covariant)
                                                             for patient in
                                                             self.matched_treatment]
            values_control_matched[(covariant, 'values')] = [getattr(patient, covariant)
                                                             for patient in
                                                             self.matched_control]

            values_treated_matched[(covariant, 'mean')] = np.mean(
                [getattr(patient, covariant) for patient in self.matched_treatment])
            values_control_matched[(covariant, 'mean')] = np.mean(
                [getattr(patient, covariant) for patient in self.matched_control])

            values_treated_matched[(covariant, 'std')] = np.std(
                [getattr(patient, covariant) for patient in self.matched_treatment])
            values_control_matched[(covariant, 'std')] = np.std(
                [getattr(patient, covariant) for patient in self.matched_control])

            values_treated_matched[(covariant, 'range')] = (np.max(
                [getattr(patient, covariant) for patient in self.matched_treatment]),
                                                            np.min(
                                                                [getattr(patient,
                                                                         covariant) for
                                                                 patient in
                                                                 self.matched_treatment]))
            values_control_matched[(covariant, 'range')] = (np.max(
                [getattr(patient, covariant) for patient in self.matched_control]),
                                                            np.min(
                                                                [getattr(patient,
                                                                         covariant) for
                                                                 patient in
                                                                 self.matched_control]))

        table_data = [["Characteristic", "Control Group", "Treatment Group"]]
        for covariant in covariates:
            table_data.append(["", "", ""])
            table_data.append([covariant, "", ""])
            for characteristic in characteristics:
                table_data.append(
                    [characteristic, values_control_matched[(covariant, characteristic)],
                     values_treated_matched[(covariant, characteristic)]])


        # Calculate Pearson correlation coefficient and p-value
        for i, j in itertools.product(range(len(covariates)), range(len(covariates))):
            if i <= j:
                continue
            table_data.append(["", "", ""])
            table_data.append([covariates[i] + ' vs ' + covariates[j], '', ''])
            correlation_coefficient_control, p_value_control = pearsonr(
                values_control_matched[(covariates[i], 'values')],
                values_control_matched[(covariates[j], 'values')])
            correlation_coefficient_treatment, p_value_treatment = pearsonr(
                values_treated_matched[(covariates[i], 'values')],
                values_treated_matched[(covariates[j], 'values')])
            table_data.append(
                ['correlation_coefficient', correlation_coefficient_control,
                 correlation_coefficient_treatment])
            table_data.append(['p_value', p_value_control, p_value_treatment])

        # Print the table
        print(tabulate(table_data, headers="firstrow", tablefmt="csv"))
