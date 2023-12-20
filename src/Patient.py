import numpy as np


class Patient:
    def __init__(self, p_vaccinated_base, age_range,
                 alpha_beta_ratio, region_params, p_infection_base, vaccination_window,
                 immunization_duration, index,
                 followup_duration):

        self.objective_day = None
        self.followup_infected = False
        self.p_infection = None
        self.index = index
        self.p_vaccinated_base = p_vaccinated_base
        self.age_range = age_range
        self.alpha_beta_ratio = alpha_beta_ratio
        self.region_params = region_params
        self.p_infection_base = p_infection_base
        self.vaccination_window = vaccination_window
        self.immunization_duration = immunization_duration
        self.followup_duration = followup_duration
        self.immunization_progress = [1] * self.followup_duration
        self.infection_days = [0] * self.followup_duration
        self.infected = False
        self.infection_time = None
        self.vaccination_days = [0] * self.followup_duration
        self.vaccinated = None
        self.vaccine_time = None
        self.ps = None
        self.matched_patient = None

        self.age = np.random.uniform(*self.age_range) / 100

        alpha = self.age * self.alpha_beta_ratio * 100
        beta = 100 - alpha
        self.past = np.random.beta(alpha, beta)

        self.region = np.random.uniform(*self.region_params)

        self.generate_trail_arm()

        self.generate_infection_progress()

        self.condition = np.random.binomial(1, self.past) if self.infected else None

    def still_at_risk(self, day):  # add whats happning if the followup ended
        return sum(self.infection_days[:day]) == 0

    def generate_trail_arm(self):
        age_past_region_effect = self.age * self.past * (1 - self.region) / 2
        self.p_vaccinated = min(max(0, self.p_vaccinated_base * age_past_region_effect), 1)
        self.vaccination_days = np.random.binomial(1, self.p_vaccinated,
                                                   size=self.followup_duration)
        self.vaccinated = (sum(self.vaccination_days) >= 1)

        if self.vaccinated:
            self.vaccine_time = next(
                (index for index, item in enumerate(self.vaccination_days) if item != 0),
                None)
            self.generate_immunization_progress()

    def generate_immunization_progress(self ):
        stages = [(0, 14), (14, 21), (21, 28), (28, self.followup_duration)]
        protection_levels = [0.05, 0.46, 0.60, 0.92]
        immunization_progress = []
        for stage, protection_level in zip(stages, protection_levels):
            immunization_progress.extend(
                1 - np.random.binomial(1, protection_level) for _ in range(*stage))

        self.immunization_progress = immunization_progress

    def generate_infection_progress(self):
        age_region_effect = self.age * self.region
        self.p_infection = self.p_infection_base * (1 + age_region_effect)

        self.infection_days = np.random.binomial(1, self.p_infection,
                                                 size=self.followup_duration)
        self.infection_days = [i * j for i, j in
                               zip(self.infection_days, self.immunization_progress)]

        self.infected = (sum(self.infection_days) >= 1)  # infected in 42 days

        if self.infected:
            self.infection_time = next(
                (index for index, item in enumerate(self.infection_days) if item != 0),
                None)

    def get_survival_data(self):
        censored = 0
        calendrical_followup_start = self.vaccine_time if self.vaccinated else self.matched_patient.vaccine_time
        calendrical_followup_end = self.followup_duration
        if self.infected:
            calendrical_followup_end, censored = min((self.followup_duration, 0), (
                calendrical_followup_start + self.infection_time, 1), key=lambda x: x[0])

        return calendrical_followup_end - calendrical_followup_start, censored, self

    def infected_while_followup(self):
        if not self.infected:
            return False

        if self.vaccinated:
            self.objective_day = self.vaccine_time + self.infection_time  # infection calanderic day
        else:
            self.objective_day = self.matched_patient.vaccine_time + self.infection_time  # infection calanderic day

        self.followup_infected = (sum(self.infection_days[:self.objective_day]) >= 1)
        return self.followup_infected

