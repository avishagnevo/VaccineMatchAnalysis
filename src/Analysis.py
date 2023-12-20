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
        self.kmf_vaccinated = None
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

        # Fit Kaplan-Meier survival curves
        self.fit_kaplan_meier()

    def fit_kaplan_meier(self):

        ci_labels = ('lower_bound', 'upper_bound')

        # Fit Kaplan-Meier survival curves
        self.kmf_vaccinated = KaplanMeierFitter()
        self.kmf_vaccinated.fit(self.event_time_treatment,
                                event_observed=self.censored_treatment,
                                ci_labels=ci_labels)

        self.kmf_unvaccinated = KaplanMeierFitter()
        self.kmf_unvaccinated.fit(self.event_time_control,
                                  event_observed=self.censored_control,
                                  ci_labels=ci_labels)

    def visualize_kaplan_meier_curves(self):
        kmf_vaccinated = self.kmf_vaccinated
        kmf_unvaccinated = self.kmf_unvaccinated

        # Plot survival curves with separate subplots
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))

        # Plot for vaccinated group
        kmf_vaccinated.plot_survival_function(ax=ax1, color='steelblue',
                                              label='Treatment', at_risk_counts=True)
        ax1.set_title('Survival Curve for Vaccinated Group')
        ax1.set_xlabel('Days')
        ax1.set_ylabel('Survival Probability')
        ax1.set_ylim(0.5, 1)
        ax1.legend()

        # Plot for unvaccinated group
        kmf_unvaccinated.plot_survival_function(ax=ax2, color='indianred',
                                                label='Control', at_risk_counts=True)
        ax2.set_title('Survival Curve for Unvaccinated Group')
        ax2.set_xlabel('Days')
        ax2.set_ylabel('Survival Probability')
        ax2.set_ylim(0.5, 1)
        ax2.legend()

        plt.show()

    def vaccine_effectivness_table(self):
        table_data = [["Stage", "Control Risks", "Treatment Risks", "Risk Ratio",
                       "Vaccine Effectiveness", "CI (95%)"]]

        stages = [(0, 14), (14, 21), (21, 28), (28, self.followup_duration)]
        for stage in stages:
            kmf_start_v = self.kmf_vaccinated.cumulative_density_at_times(stage[0]).loc[
                stage[0]]
            kmf_end_v = self.kmf_vaccinated.cumulative_density_at_times(stage[1]).loc[
                stage[1]]

            kmf_start_u = self.kmf_unvaccinated.cumulative_density_at_times(stage[0]).loc[
                stage[0]]
            kmf_end_u = self.kmf_unvaccinated.cumulative_density_at_times(stage[1]).loc[
                stage[1]]

            stage_risk_ratio = (kmf_end_v - kmf_start_v) / (kmf_end_u - kmf_start_u)
            rr_ci = self.get_risk_ratio_ci(stage)

            table_data.append(
                [stage, f"{kmf_end_u:.03} -> {kmf_start_u:.03}",
                 f"{kmf_end_v:.03} -> {kmf_start_v:.03}",
                 f"{stage_risk_ratio:.03}", f"{1 - stage_risk_ratio:.03}", rr_ci])
        # Print the table
        print(
            'Estimated Vaccine Effectivness against Outcome during Four Time Periods:')
        print(tabulate(table_data, headers="firstrow", tablefmt="csv"))
        print()

    def get_risk_ratio_ci(self, stage):
        kmf_ci_start_v = np.array(
            [self.kmf_vaccinated.confidence_interval_['lower_bound'][stage[0]],
             self.kmf_vaccinated.confidence_interval_['upper_bound'][stage[0]]])
        kmf_ci_end_v = np.array(
            [self.kmf_vaccinated.confidence_interval_['lower_bound'][stage[1]],
             self.kmf_vaccinated.confidence_interval_['upper_bound'][stage[1]]])

        kmf_ci_start_u = np.array(
            [self.kmf_unvaccinated.confidence_interval_['lower_bound'][stage[0]],
             self.kmf_unvaccinated.confidence_interval_['upper_bound'][stage[0]]])
        kmf_ci_end_u = np.array(
            [self.kmf_unvaccinated.confidence_interval_['lower_bound'][stage[1]],
             self.kmf_unvaccinated.confidence_interval_['upper_bound'][stage[1]]])

        kmf_log_ci_v = np.log(kmf_ci_start_v - kmf_ci_end_v)
        kmf_log_ci_u = np.log(kmf_ci_start_u - kmf_ci_end_u)

        rr_log_ci = kmf_log_ci_v - kmf_log_ci_u

        return ["{:.4f}".format(1 - bound) for bound in np.exp(rr_log_ci)]

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

            values_treated_matched[(covariant,
                                    'mean')] = f"{np.mean([getattr(patient, covariant) for patient in self.matched_treatment]):.03}"
            values_control_matched[(covariant,
                                    'mean')] = f"{np.mean([getattr(patient, covariant) for patient in self.matched_control]):.03}"

            values_treated_matched[(covariant,
                                    'std')] = f"{np.std([getattr(patient, covariant) for patient in self.matched_treatment]):.03}"
            values_control_matched[(covariant,
                                    'std')] = f"{np.std([getattr(patient, covariant) for patient in self.matched_control]):.03}"

            values_treated_matched[(covariant, 'range')] = (
                f"{np.min([getattr(patient, covariant) for patient in self.matched_treatment]):.03}",
                f"{np.max([getattr(patient, covariant) for patient in self.matched_treatment]):.03}")
            values_control_matched[(covariant, 'range')] = (
                f"{np.min([getattr(patient, covariant) for patient in self.matched_control]):.03}",
                f"{np.max([getattr(patient, covariant) for patient in self.matched_control]):.03}")

        table_data = [
            ["Characteristic", f"Control Group (N= {len(self.matched_control)})",
             f"Treatment Group (N= {len(self.matched_treatment)})"]]
        for covariant in covariates:
            table_data.append([covariant, "", ""])
            table_data.extend(
                [characteristic, values_control_matched[(covariant, characteristic)],
                 values_treated_matched[(covariant, characteristic)]] for characteristic
                in characteristics)

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
                ['correlation_coefficient', f"{correlation_coefficient_control:.03}",
                 f"{correlation_coefficient_treatment:.03}"])

        # Print the table
        print('Table 1:')
        print(tabulate(table_data, headers="firstrow", tablefmt="csv"))
        print()

    def visualize_vaccine_effectivness(self):
        days = range(self.followup_duration)
        kmf_v = self.kmf_vaccinated.cumulative_density_at_times(days)
        kmf_u = self.kmf_unvaccinated.cumulative_density_at_times(days)
        daily_risk_ratios = []
        daily_vaccine_effectivness = []
        for time in days[:-1]:
            kmf_start_v = kmf_v.loc[time]
            kmf_end_v = kmf_v.loc[time + 1]
            kmf_start_u = kmf_u.loc[time]
            kmf_end_u = kmf_u.loc[time + 1]
            RR = (kmf_end_v - kmf_start_v) / (kmf_end_u - kmf_start_u)
            daily_risk_ratios.append(RR)
            daily_vaccine_effectivness.append(1 - RR)

        plt.figure(figsize=(12, 6))
        plt.plot(days[:-1], daily_vaccine_effectivness, marker='p',
                 c='black', label='Vaccine Effectivness')
        plt.title('Vaccine Effectiveness Over Time')
        plt.xlabel('Time')
        plt.ylabel('1 - RR')
        plt.legend()
        plt.text(0.5, -0.1, '* Missing data points is for times that the RR is '
                            'undefined : infinitesimal to no change at control group / both groups together',
                 ha='center', va='center',
                 transform=plt.gca().transAxes)
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
