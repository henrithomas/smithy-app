// Chart.js defaults
// Chart.defaults.font.size = 15
// Chart.defaults.plugins.title.font.size = 20
// Chart.defaults.plugins.legend.position = 'bottom'

var times = [];
var costs = [];

// SETUP
const gibson_bool = !!document.getElementById('gibson-dashboard');
const goldengate_bool = !!document.getElementById('goldengate-dashboard')
const slic_bool = !!document.getElementById('slic-dashboard')
const pcr_bool = !!document.getElementById('pcr-dashboard')
const biobricks_bool = !!document.getElementById('biobricks-dashboard')
// data
// Gibson
if (gibson_bool) {
    const gibson_time_data = JSON.parse(JSON.parse(document.getElementById("gibson_times").textContent));
    const gibson_cost_data = JSON.parse(JSON.parse(document.getElementById("gibson_costs").textContent));
    const gibson_risk_data = JSON.parse(JSON.parse(document.getElementById("gibson_risks").textContent)); 

    times.push(gibson_time_data.total);
    costs.push(gibson_cost_data.total);
}

// Golden Gate
if (goldengate_bool) {
    const goldengate_time_data = JSON.parse(JSON.parse(document.getElementById("goldengate_times").textContent));
    const goldengate_cost_data = JSON.parse(JSON.parse(document.getElementById("goldengate_costs").textContent));
    const goldengate_risk_data = JSON.parse(JSON.parse(document.getElementById("goldengate_risks").textContent));

    times.push(goldengate_time_data.total);
    costs.push(goldengate_cost_data.total);
}

// SLIC
if (slic_bool) {
    const slic_time_data = JSON.parse(JSON.parse(document.getElementById("slic_times").textContent));
    const slic_cost_data = JSON.parse(JSON.parse(document.getElementById("slic_costs").textContent));
    const slic_risk_data = JSON.parse(JSON.parse(document.getElementById("slic_risks").textContent));

    times.push(slic_time_data.total);
    costs.push(slic_cost_data.total);
}

// PCR
if (pcr_bool) {
    const pcr_time_data = JSON.parse(JSON.parse(document.getElementById("pcr_times").textContent));
    const pcr_cost_data = JSON.parse(JSON.parse(document.getElementById("pcr_costs").textContent));
    const pcr_risk_data = JSON.parse(JSON.parse(document.getElementById("pcr_risks").textContent));

    times.push(pcr_time_data.total);
    costs.push(pcr_cost_data.total);
}

// BioBricks
if (biobricks_bool) {
    const biobricks_time_data = JSON.parse(JSON.parse(document.getElementById("biobricks_times").textContent));
    const biobricks_cost_data = JSON.parse(JSON.parse(document.getElementById("biobricks_costs").textContent));
    const biobricks_risk_data = JSON.parse(JSON.parse(document.getElementById("biobricks_risks").textContent));

    times.push(biobricks_time_data.total);
    costs.push(biobricks_cost_data.total);
}

// find max values here for time and cost

// alter datasets 

if (gibson_bool) {
    gibson_time_chart(gibson_time_data.times, gibson_time_data.types, gibson_time_data.total, "gibson-time", "gibson-time-sum");
    gibson_cost_chart(gibson_cost_data.costs, gibson_cost_data.types, gibson_cost_data.total, "gibson-cost");
    gibson_risk_chart(gibson_risk_data.risks, gibson_risk_data.types, "gibson-risk");
}

if (goldengate_bool) {
    goldengate_time_chart(goldengate_time_data.times, goldengate_time_data.types, goldengate_time_data.total, "goldengate-time", "goldengate-time-sum");
    goldengate_cost_chart(goldengate_cost_data.costs, goldengate_cost_data.types, goldengate_cost_data.total, "goldengate-cost");
    goldengate_risk_chart(goldengate_risk_data.risks, goldengate_risk_data.types, "goldengate-risk");
}

if (slic_bool) {
    slic_time_chart(slic_time_data.times, slic_time_data.types, slic_time_data.total, "slic-time", "slic-time-sum");
    slic_cost_chart(slic_cost_data.costs, slic_cost_data.types, slic_cost_data.total, "slic-cost");
    slic_risk_chart(slic_risk_data.risks, slic_risk_data.types, "slic-risk");
}

if (pcr_bool) {
    pcr_time_chart(pcr_time_data.times, pcr_time_data.types, pcr_time_data.total, "pcr-time", "pcr-time-sum");
    pcr_cost_chart(pcr_cost_data.costs, pcr_cost_data.types, pcr_cost_data.total, "pcr-cost");
    pcr_risk_chart(pcr_risk_data.risks, pcr_risk_data.types, "pcr-risk");
}

if (biobricks_bool) {
    biobricks_time_chart(biobricks_time_data.times, biobricks_time_data.types, biobricks_time_data.total, "biobricks-time", "biobricks-time-sum");
    biobricks_cost_chart(biobricks_cost_data.costs, biobricks_cost_data.types, biobricks_cost_data.total, "biobricks-cost");
    biobricks_risk_chart(biobricks_risk_data.risks, biobricks_risk_data.types, "biobricks-risk");
}