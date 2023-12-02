import matplotlib.pyplot as plt
from Patient import Patient
import numpy as np
from sklearn.linear_model import LogisticRegression


class Experiment:
    def __init__(self, num_samples, p_vaccinated_base, age_range,
                 alpha_beta_ratio, region_params, p_infection_base, vaccination_window,
                 immunization_duration, followup_duration):

        self.num_samples = num_samples
        self.p_vaccinated_base = p_vaccinated_base
        self.age_range = age_range
        self.alpha_beta_ratio = alpha_beta_ratio
        self.region_params = region_params
        self.p_infection_base = p_infection_base
        self.vaccination_window = vaccination_window
        self.immunization_duration = immunization_duration
        self.followup_duration = followup_duration
        self.data = []
        self.vaccination_values = []
        self.age_values = []
        self.past_values = []
        self.region_values = []
        self.infection_prob_values = []
        self.infection_values = []

        # For analysis
        self.vaccinated_infected = []
        self.vaccinated_not_infected = []
        self.not_vaccinated_infected = []
        self.not_vaccinated_not_infected = []
        self.infected_count = None
        self.vaccinated_count = None

        self.patients = []
        self.propensity_scores = None

    def generate_data(self):
        for i in range(self.num_samples):

            patient = Patient(self.p_vaccinated_base, self.age_range,
                              self.alpha_beta_ratio, self.region_params,
                              self.p_infection_base, self.vaccination_window,
                              self.immunization_duration, i, self.followup_duration)

            self.age_values.append(patient.age)
            self.past_values.append(patient.past)
            self.region_values.append(patient.region)
            self.vaccination_values.append(patient.vaccinated)
            self.infection_prob_values.append(patient.p_infection)
            self.infection_values.append(patient.infection_time is not None)

            if patient.vaccinated and patient.infected:
                self.vaccinated_infected.append(patient)
            elif patient.vaccinated:
                self.vaccinated_not_infected.append(patient)
            elif patient.infected:
                self.not_vaccinated_infected.append(patient)
            else:
                self.not_vaccinated_not_infected.append(patient)

            self.data.append(
                [patient.age, patient.past, patient.region, patient.vaccinated,
                 patient.vaccine_time, patient.infected, patient.infection_time,
                 patient.condition])

            self.patients.append(patient)

        self.vaccinated_count = sum(self.vaccination_values)
        self.infected_count = sum(self.infection_values)

    def set_propensity_scores(self):
        # Build a logistic regression model
        X = np.array([self.age_values, self.past_values, self.region_values]).T
        y = np.array(self.vaccination_values)
        model = LogisticRegression()
        model.fit(X, y)

        # Get propensity scores
        self.propensity_scores = model.predict_proba(X)[:, 1]

        for patient, ps in zip(self.patients, self.propensity_scores):
            patient.ps = ps


    def visualize_data(self):
        # Visualizations
        plt.figure(figsize=(12, 12))

        # Age Distribution
        plt.subplot(2, 3, 1)
        plt.hist(self.age_values, bins=20, edgecolor='k')
        plt.xlabel('Age')
        plt.ylabel('Frequency')
        plt.title('Age Distribution')

        # Past Distribution
        plt.subplot(2, 3, 2)
        plt.hist(self.past_values, bins=20, edgecolor='k')
        plt.xlabel('Past Risk')
        plt.ylabel('Frequency')
        plt.title('Past Distribution')

        # Region Distribution
        plt.subplot(2, 3, 3)
        plt.hist(self.region_values, bins=20, edgecolor='k')
        plt.xlabel('Region Risk')
        plt.ylabel('Frequency')
        plt.title('Region Distribution')

        # Age-Past Scatter Plot
        plt.subplot(2, 3, 4)
        plt.scatter(self.age_values, self.past_values, marker='o', alpha=0.5)
        plt.xlabel('Age')
        plt.ylabel('Past Risk')
        plt.title('Age vs. Past Risk')

        # Past-infected Scatter Plot
        plt.subplot(2, 3, 5)
        plt.scatter(self.past_values, self.infection_prob_values, marker='o', alpha=0.5)
        plt.xlabel('Past')
        plt.ylabel('Chance of infection')
        plt.title('Past vs Infected probability values')

        # Region-infected Scatter Plot
        plt.subplot(2, 3, 6)
        plt.scatter(self.past_values, self.infection_prob_values, marker='o', alpha=0.5)
        plt.xlabel('Region')
        plt.ylabel('Chance of infection')
        plt.title('Region vs Infected probability values')

        # Region-Past-Infection probability Scatter Plot
        plt.figure(figsize=(12, 6))
        plt.scatter(self.region_values, self.past_values, c=self.infection_prob_values,
                    cmap='coolwarm',
                    alpha=0.5)
        plt.xlabel('Region Risk')
        plt.ylabel('Past Risk')
        plt.title('Effect of Region and Past Risk on Daily Infection Probability')
        cbar = plt.colorbar()
        cbar.set_label('Daily Infection Probability')

        # Pie chart: Vaccinated and Infected
        plt.figure(figsize=(12, 6))
        labels = 'Treatment Infected', 'Treatment Not Infected', 'Control Infected', 'Control Not Infected'
        sizes = [len(self.vaccinated_infected), len(self.vaccinated_not_infected),
                 len(self.not_vaccinated_infected),
                 len(self.not_vaccinated_not_infected)]
        colors = ['steelblue', 'lightskyblue', 'indianred', 'lightcoral']
        plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=140)
        plt.axis('equal')
        plt.title('Vaccine and Infection chart')

        plt.tight_layout()
        plt.show()

        ''' 
        More Visualizations:
        # Region-Past-Infection Scatter Plot
        plt.figure(figsize=(12, 6))
        plt.scatter(self.region_values, self.past_values, c=self.infection_values,
                    cmap='bwr', alpha=0.5)
        plt.xlabel('Region Risk')
        plt.ylabel('Past Risk')
        plt.title('Effect of Region and Past Risk on Infection')
        infected_patch = mpatches.Patch(color='red', label='Infected')
        not_infected_patch = mpatches.Patch(color='blue', label='Not Infected')
        plt.legend(handles=[infected_patch, not_infected_patch])


        # Pie chart: Vaccinated vs. Not Vaccinated
        vaccinated_labels = 'Vaccinated', 'Not Vaccinated'
        vaccinated_sizes = [self.vaccinated_count,
                            self.num_samples - self.vaccinated_count]
        vaccinated_colors = ['steelblue', 'indianred']
        plt.figure(figsize=(8, 4))
        plt.subplot(1, 2, 1)
        plt.pie(vaccinated_sizes, labels=vaccinated_labels, colors=vaccinated_colors,
                autopct='%1.1f%%', startangle=140)
        plt.axis('equal')
        plt.title('Vaccinated vs. Not Vaccinated')

        # Pie chart: Infected vs. Not Infected
        infected_labels = 'Infected', 'Not Infected'
        infected_sizes = [self.infected_count, self.num_samples - self.infected_count]
        infected_colors = ['crimson', 'cadetblue']
        plt.subplot(1, 2, 2)
        plt.pie(infected_sizes, labels=infected_labels, colors=infected_colors,
                autopct='%1.1f%%', startangle=140)
        plt.axis('equal')
        plt.title('Infected vs. Not Infected')
        '''

    def visualize_propensity_scores(self):
        # propensity score visualization

        # Region-Past-propensity score Scatter Plot
        plt.figure(figsize=(12, 6))
        plt.subplot(1, 3, 1)
        plt.scatter(self.region_values, self.past_values, c=self.propensity_scores,
                    cmap='coolwarm',
                    alpha=0.5)
        plt.xlabel('Region Risk')
        plt.ylabel('Past Risk')
        plt.title('Effect of Region and Past Risk on PS')

        plt.subplot(1, 3, 2)
        plt.scatter(self.region_values, self.age_values, c=self.propensity_scores,
                    cmap='coolwarm',
                    alpha=0.5)
        plt.xlabel('Region Risk')
        plt.ylabel('Age')
        plt.title('Effect of Region and Age on PS')

        plt.subplot(1, 3, 3)
        plt.scatter(self.age_values, self.past_values, c=self.propensity_scores,
                    cmap='coolwarm',
                    alpha=0.5)
        plt.xlabel('Age Risk')
        plt.ylabel('Past Risk')
        plt.title('Effect of Age and Past Risk on PS')
        cbar = plt.colorbar()
        cbar.set_label('Propensity Score')

        plt.tight_layout()
        plt.show()
