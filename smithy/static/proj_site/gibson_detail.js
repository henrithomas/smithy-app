const time_data = JSON.parse(JSON.parse(document.getElementById("times").textContent));
const cost_data = JSON.parse(JSON.parse(document.getElementById("costs").textContent));
const risk_data = JSON.parse(JSON.parse(document.getElementById("risks").textContent));

gibson_time_chart(time_data.times, time_data.types, time_data.total, "gibson-time", "gibson-time-sum");
gibson_cost_chart(cost_data.costs, cost_data.types, cost_data.total, "gibson-cost");
gibson_risk_chart(risk_data.risks, risk_data.types, "gibson-risk");