const time_data = JSON.parse(JSON.parse(document.getElementById("times").textContent));
const cost_data = JSON.parse(JSON.parse(document.getElementById("costs").textContent));
const risk_data = JSON.parse(JSON.parse(document.getElementById("risks").textContent));

pcr_time_chart(time_data.times, time_data.types, time_data.total, "pcr-time", "pcr-time-sum");
pcr_cost_chart(cost_data.costs, cost_data.types, cost_data.total, "pcr-cost");
pcr_risk_chart(risk_data.risks, risk_data.types, "pcr-risk");