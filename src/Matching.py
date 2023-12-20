from Experiment import *


class Matching:
    def __init__(self, experiment, caliper=0.1):
        self.patients = experiment.patients
        self.caliper = caliper
        self.treatment = set()
        self.control = set()
        self.vaccination_window = experiment.vaccination_window
        self.matched_treatment = None
        self.matched_control = None
        self.day = 0
        self.experiment = experiment

        for patient in self.patients:
            if patient.vaccination_days[0]:
                self.treatment.add(patient)
            else:
                self.control.add(patient)

    def update_treatment_control(self, day):
    # sequential continues updating of treatment control sets
        control = self.control.copy()
        for patient in control:
            if patient.vaccination_days[day]:  # in vaccination period
                self.control.remove(patient)
                self.treatment.add(patient)
                if patient.matched_patient is not None:
                    self.treatment.remove(patient.matched_patient)
                    patient.matched_patient = None

    def run_daily_matching(self, day):
        self.day = day
        if day < self.experiment.followup_duration:
            self.update_treatment_control(day)

        # Match the treatment (exposed group) to control using a narrow caliper
        for t_patient in self.treatment:
            if t_patient.matched_patient is not None:  # if a patient has a match
                continue

            followup_day = self.day - t_patient.vaccine_time
            if not t_patient.still_at_risk(
                    day=followup_day):  # doing matching only on the healthy
                continue

            treated_propensity = t_patient.ps
            best_control_index = None
            best_distance = np.inf

            for c_patient in self.control:
                if not c_patient.still_at_risk(
                        day=followup_day):  # doing matching only on the healthy
                    continue
                control_index = c_patient.index
                if c_patient.matched_patient is None:  # Ensure control individuals are not used more than once
                    control_propensity = c_patient.ps
                    distance = abs(treated_propensity - control_propensity)

                    if distance < self.caliper and distance < best_distance:
                        best_control_index = control_index
                        best_distance = distance

            if best_control_index is not None:
                c_patient = self.patients[best_control_index]
                t_patient.matched_patient = c_patient
                c_patient.matched_patient = t_patient

    def run_matching(self):
        # day is day of follow up
        for day in range(
                self.vaccination_window):  # longer? what about a patient that was vaccinated at 31.1.2020
            matched_treatment = []
            matched_control = []

            self.run_daily_matching(day)

            for patient in self.treatment:  # only patients survived in trail
                if patient.matched_patient is not None:
                    matched_treatment.append(patient)
                    matched_control.append(patient.matched_patient)

            self.matched_treatment = matched_treatment
            self.matched_control = matched_control

    def visualize_matching(self):
        # Visualize the matching and propensity scores
        treatment_color = ['steelblue' for _ in self.matched_treatment]
        control_color = ['indianred' for _ in self.matched_control]
        ps_treatment = [patient.ps for patient in self.matched_treatment]
        ps_control = [patient.ps for patient in self.matched_control]
        length = len(self.matched_treatment)

        plt.figure(figsize=(12, 6))
        plt.scatter(range(length), ps_treatment, c=treatment_color, marker='_', label='Treatment')
        plt.scatter(range(length), ps_control, c=control_color, marker= '|', label='Control')
        plt.xlabel('Individuals')
        plt.ylabel('Propensity Score')
        plt.title('Matching and Propensity Scores')
        plt.legend()
        plt.show()

        # Create box plots for 0.1 propensity score intervals
        propensity_intervals = np.arange(0, 1.1, 0.1)
        interval_distances = [[] for _ in range(len(propensity_intervals) - 1)]

        for i in range(len(self.matched_treatment)):
            treated_patient = self.matched_treatment[i]
            control_patient = treated_patient.matched_patient
            propensity_diff = abs(treated_patient.ps - control_patient.ps)

            for j in range(len(propensity_intervals) - 1):
                if propensity_intervals[j] <= treated_patient.ps < propensity_intervals[
                    j + 1]:
                    interval_distances[j].append(propensity_diff)

            # Plot box plots
        plt.figure(figsize=(12, 6))
        boxplot = plt.boxplot(interval_distances, labels=[
            f"{round(propensity_intervals[i], 1)}-{round(propensity_intervals[i + 1], 1)}"
            for i in range(len(propensity_intervals) - 1)],
                    meanline=True, showmeans=True, meanprops=dict(color='blue'))

        plt.xlabel('Propensity Score')
        plt.ylabel('Distance Between Matched Patients')
        plt.title(
            'Distance Distribution Between Matched Patients by Propensity Score Intervals')
        plt.legend([boxplot["means"][0]], ['Mean'], loc='upper right')
        plt.show()

    def visualize_means_diff(self):
        covariates = ['age', 'past', 'region']
        means_treated_matched = {}
        means_control_matched = {}
        covariant_values_matched = {}
        std_diff_after = {}
        means_treated = {}
        means_control = {}
        covariant_values = {}
        std_diff_before = {}

        for covariant in covariates:
            means_treated_matched[covariant] = np.mean(
                [getattr(patient, covariant) for patient in self.matched_treatment])
            means_control_matched[covariant] = np.mean(
                [getattr(patient, covariant) for patient in self.matched_control])
            covariant_values_matched[covariant] = [getattr(patient, covariant)
                                                   for patient in
                                                   self.matched_treatment + self.matched_control]

            std_diff_after[covariant] = abs(
                means_treated_matched[f'{covariant}'] - means_control_matched[
                    f'{covariant}'])

            means_treated[covariant] = np.mean(
                [getattr(patient, covariant) for patient in self.treatment])
            means_control[covariant] = np.mean(
                [getattr(patient, covariant) for patient in self.control])
            covariant_values[covariant] = [getattr(patient, covariant)
                                           for patient in
                                           self.matched_treatment + self.matched_control]

            std_diff_before[covariant] = abs(
                means_treated[f'{covariant}'] - means_control[f'{covariant}'])

        # Plot standardized differences in means before and after matching
        std_diffs_before = [std_diff_before[covariant] for covariant in covariates]
        std_diffs_after = [std_diff_after[covariant] for covariant in covariates]



        '''
        plt.figure(figsize=(10, 6))
        plt.bar(covariates, std_diffs_before, width=0.1, label='Before Matching',
                align='center')
        plt.bar(covariates, std_diffs_after, width=0.1, label='After Matching',
                align='edge', alpha=0.7)

        plt.xlabel('Covariate')
        plt.ylabel('Standardized Difference')
        plt.title('Standardized Differences in Means Before and After Matching')
        plt.legend()
        plt.show()
        '''

        plt.figure(figsize=(12, 6))
        # Plotting lollipops for Before Matching
        plt.stem(covariates, std_diffs_before, markerfmt='bo', linefmt='b-', basefmt=' ')
        # Plotting lollipops for After Matching
        plt.stem(covariates, std_diffs_after, markerfmt='ro', linefmt='r-', basefmt=' ')

        plt.xlabel('Covariate')
        plt.ylabel('Standardized Mean Difference')
        plt.title('Standardized Differences in Means Before and After Matching')
        plt.legend(['Before Matching', 'After Matching'])
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.show()
