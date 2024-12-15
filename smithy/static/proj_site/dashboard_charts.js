// Chart.js defaults
// Chart.defaults.font.size = 15
// Chart.defaults.plugins.title.font.size = 20
// Chart.defaults.plugins.legend.position = 'bottom'

// data: {
//     costs: [x, y, z,...],
//     total: a,
//     types: [a, b, c,...]
// }

function offset_check(total, max) {
    if (total < max) {
        return [true, max - total];
    }
    else {
        return [false, 0]
    }
}

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
    var gibson_time_data = JSON.parse(JSON.parse(document.getElementById("gibson_times").textContent));
    var gibson_cost_data = JSON.parse(JSON.parse(document.getElementById("gibson_costs").textContent));
    var gibson_risk_data = JSON.parse(JSON.parse(document.getElementById("gibson_risks").textContent)); 

    times.push(gibson_time_data.total);
    costs.push(gibson_cost_data.total);
}

// Golden Gate
if (goldengate_bool) {
    var goldengate_time_data = JSON.parse(JSON.parse(document.getElementById("goldengate_times").textContent));
    var goldengate_cost_data = JSON.parse(JSON.parse(document.getElementById("goldengate_costs").textContent));
    var goldengate_risk_data = JSON.parse(JSON.parse(document.getElementById("goldengate_risks").textContent));

    times.push(goldengate_time_data.total);
    costs.push(goldengate_cost_data.total);
}

// SLIC
if (slic_bool) {
    var slic_time_data = JSON.parse(JSON.parse(document.getElementById("slic_times").textContent));
    var slic_cost_data = JSON.parse(JSON.parse(document.getElementById("slic_costs").textContent));
    var slic_risk_data = JSON.parse(JSON.parse(document.getElementById("slic_risks").textContent));

    times.push(slic_time_data.total);
    costs.push(slic_cost_data.total);
}

// PCR
if (pcr_bool) {
    var pcr_time_data = JSON.parse(JSON.parse(document.getElementById("pcr_times").textContent));
    var pcr_cost_data = JSON.parse(JSON.parse(document.getElementById("pcr_costs").textContent));
    var pcr_risk_data = JSON.parse(JSON.parse(document.getElementById("pcr_risks").textContent));

    times.push(pcr_time_data.total);
    costs.push(pcr_cost_data.total);
}

// BioBricks
if (biobricks_bool) {
    var biobricks_time_data = JSON.parse(JSON.parse(document.getElementById("biobricks_times").textContent));
    var biobricks_cost_data = JSON.parse(JSON.parse(document.getElementById("biobricks_costs").textContent));
    var biobricks_risk_data = JSON.parse(JSON.parse(document.getElementById("biobricks_risks").textContent));

    times.push(biobricks_time_data.total);
    costs.push(biobricks_cost_data.total);
}

// find max values here for time and cost
var timeMax = Math.ceil(Math.max.apply(null, times));
var costMax = Math.max.apply(null, costs);


if (gibson_bool) {
    let [gibson_offset, gibson_offset_amt] = offset_check(gibson_cost_data.total, costMax);

    gibson_time_chart(gibson_time_data.times, gibson_time_data.types, gibson_time_data.total, "gibson-time", "gibson-time-sum", timeMax);
    gibson_cost_chart(gibson_cost_data.costs, gibson_cost_data.types, gibson_cost_data.total, "gibson-cost", gibson_offset, gibson_offset_amt);
    gibson_risk_chart(gibson_risk_data.risks, gibson_risk_data.types, "gibson-risk");
}

if (goldengate_bool) {
    let [goldengate_offset, goldengate_offset_amt] = offset_check(goldengate_cost_data.total, costMax);

    goldengate_time_chart(goldengate_time_data.times, goldengate_time_data.types, goldengate_time_data.total, "goldengate-time", "goldengate-time-sum", timeMax);
    goldengate_cost_chart(goldengate_cost_data.costs, goldengate_cost_data.types, goldengate_cost_data.total, "goldengate-cost", goldengate_offset, goldengate_offset_amt);
    goldengate_risk_chart(goldengate_risk_data.risks, goldengate_risk_data.types, "goldengate-risk");
}

if (slic_bool) {
    let [slic_offset, slic_offset_amt] = offset_check(slic_cost_data.total, costMax);

    slic_time_chart(slic_time_data.times, slic_time_data.types, slic_time_data.total, "slic-time", "slic-time-sum", timeMax);
    slic_cost_chart(slic_cost_data.costs, slic_cost_data.types, slic_cost_data.total, "slic-cost", slic_offset, slic_offset_amt);
    slic_risk_chart(slic_risk_data.risks, slic_risk_data.types, "slic-risk");
}

if (pcr_bool) {
    let [pcr_offset, pcr_offset_amt] = offset_check(pcr_cost_data.total, costMax);

    pcr_time_chart(pcr_time_data.times, pcr_time_data.types, pcr_time_data.total, "pcr-time", "pcr-time-sum", timeMax);
    pcr_cost_chart(pcr_cost_data.costs, pcr_cost_data.types, pcr_cost_data.total, "pcr-cost", pcr_offset, pcr_offset_amt);
    pcr_risk_chart(pcr_risk_data.risks, pcr_risk_data.types, "pcr-risk");
}

if (biobricks_bool) {
    let [biobricks_offset, biobricks_offset_amt] = offset_check(biobricks_cost_data.total, costMax);

    biobricks_time_chart(biobricks_time_data.times, biobricks_time_data.types, biobricks_time_data.total, "biobricks-time", "biobricks-time-sum", timeMax);
    biobricks_cost_chart(biobricks_cost_data.costs, biobricks_cost_data.types, biobricks_cost_data.total, "biobricks-cost", biobricks_offset, biobricks_offset_amt);
    biobricks_risk_chart(biobricks_risk_data.risks, biobricks_risk_data.types, "biobricks-risk");
}