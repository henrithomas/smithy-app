// Chart.js defaults
// Chart.defaults.font.size = 15
// Chart.defaults.plugins.title.font.size = 20
// Chart.defaults.plugins.legend.position = 'bottom'

// checker = !!document.getElementById("gibson-time");
// test_html = document.getElementById("testtext");
// test_html.innerHTML = "TEST";

// SETUP
// data
const gibson_times = [27, 15, 12];
const goldengate_times = [30, 42, 9];
const pcr_times = [30, 15, 26];

const gibson_costs = [65, 59, 90, 81, 56, 55];
const goldengate_costs = [65, 59, 90, 81, 161];

const data2 = [65, 59, 90, 81, 56, 55];
let cost_compare = [];

const gibson_risks = [-.9, .8, 0.1, 0.5];
const goldengate_risks = [-.5, .6, 0.3, -0.75];

const time_sum_colors = [];
const cost_sum_colors = [];
const risk_sum_colors = [];

//defaults
Chart.defaults.scales.linear.max = 81;


// const <data_setup> = {
//     labels: <labels>,
//     datasets: [
//         <data>
//     ]
// }

gibson_time_chart(gibson_times, "gibson-time", "gibson-time-sum");
gibson_cost_chart(gibson_costs, "gibson-cost");
gibson_risk_chart(gibson_risks, "gibson-risk");

goldengate_time_chart(goldengate_times, "goldengate-time", "goldengate-time-sum");
goldengate_cost_chart(goldengate_costs, "goldengate-cost");
goldengate_risk_chart(goldengate_risks, "goldengate-risk");

pcr_time_chart(pcr_times, "pcr-time", "pcr-time-sum");


// CONFIG
// const config = {
//     type: 'line',
//     data: data,
//     options: {}
//   };