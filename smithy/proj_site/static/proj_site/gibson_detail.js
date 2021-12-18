const gibson_times = [34, 6, 6];
const gibson_costs = [65, 230, 30, 24, 56, 55];
const gibson_risks = [-.9, .8, 0.1, 0.5];

gibson_time_chart(gibson_times, "gibson-time", "gibson-time-sum");
gibson_cost_chart(gibson_costs, "gibson-cost");
gibson_risk_chart(gibson_risks, "gibson-risk");