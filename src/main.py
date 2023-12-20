from Analysis import Analysis
from Experiment import Experiment
from Matching import Matching


def experiment(num_samples=50000,
               p_vaccinated_base=0.7,
               age_lower=0,
               age_upper=90,
               alpha_beta_ratio=0.5,
               region_lower=0,
               region_upper=1,
               p_infection_base=0.005,
               vaccination_window=42,
               immunization_duration=0,
               followup_duration=42
               ):

    return Experiment(num_samples, p_vaccinated_base, (age_lower, age_upper),
                      alpha_beta_ratio, [region_lower, region_upper], p_infection_base,
                      vaccination_window, immunization_duration, followup_duration)


if __name__ == '__main__':
    experiment = experiment()
    experiment.generate_data()
    experiment.visualize_data()
    experiment.set_propensity_scores()
    experiment.visualize_propensity_scores()
    matcher = Matching(experiment)
    matcher.run_matching()
    matcher.visualize_matching()
    matcher.visualize_means_diff()
    analysis = Analysis(experiment, matcher)
    analysis.create_table_1()
    analysis.visualize_kaplan_meier_curves()
    analysis.vaccine_effectivness_table()
    analysis.visualize_vaccine_effectivness()
    #analysis.cox_proportional_hazard()

