from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test

from matplotlib import pyplot as plt

from lifelines import CoxPHFitter
import numpy as np

def calculateHR(df,referenceGroupLabel,expGroupLabel,event, time):
    temp_df=df.copy()
    temp_df['group']=temp_df['group'].replace(referenceGroupLabel,0).replace(expGroupLabel,1)
    cph = CoxPHFitter().fit(temp_df, duration_col=time, event_col=event)
    # calculate p-value
    p_value = cph._compute_p_values()[0]
    # HR
    hazard_ratio=cph.hazard_ratios_.iloc[0]
    # calculate confident interval
    CI=cph.confidence_intervals_
    lower_bound=np.exp(CI.loc['group','95% lower-bound'])
    upper_bound=np.exp(CI.loc['group','95% upper-bound'])
    print("Hazard Ratio: {:.2f}".format(hazard_ratio))
    print("Confidence Interval: [{:0.2f}, {:0.2f}]".format(lower_bound, upper_bound))
    print("p-value: {:.2e}".format(p_value))
    return [p_value,hazard_ratio,lower_bound,upper_bound]



def plot(df, event, time, controlLabel="control",controlColor='r',expLabel="miR-137",expColor='g'):

    ax = plt.subplot(111)

    kmf_control = KaplanMeierFitter()
    
    groups_control = df[df.group==controlLabel]
    kmf_control.fit(groups_control[time], groups_control[event], label=controlLabel)
    ax = kmf_control.plot_survival_function(color=controlColor,show_censors=True, censor_styles={'ms': 6, 'marker': 'x'})
    control_median_survival_time=kmf_control.median_survival_time_
    print ("control median survival time:"+str(control_median_survival_time))

    kmf_exp = KaplanMeierFitter()

    groups_exp = df[df.group==expLabel]
    kmf_exp.fit(groups_exp[time], groups_exp[event], label=expLabel)
    ax = kmf_exp.plot_survival_function(color=expColor,show_censors=True, censor_styles={'ms': 6, 'marker': 'x'})
    experimental_median_survival_time=kmf_exp.median_survival_time_
    print ("exp median survival time:"+str(experimental_median_survival_time))

    add_at_risk_counts(kmf_exp, kmf_control, ax=ax)
    plt.plot()
    # Add dashed lines for median survival times
    plt.axvline(control_median_survival_time, color='grey', linestyle='--', ymax=0.5, label=f'Control Median: {control_median_survival_time:.2f}')
    plt.axvline(experimental_median_survival_time, color='grey', linestyle='--', ymax=0.5,label=f'Experimental Median: {experimental_median_survival_time:.2f}')
    

    results = logrank_test(groups_exp[time], groups_control[time], event_observed_A=groups_exp[event], event_observed_B=groups_control[event])
    #results.print_summary()
    print("log-rank test p-value: {:.2e}".format(results.p_value))