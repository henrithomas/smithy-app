function goldengate_time_chart(data, labels, sum, chart_elem, sum_elem, yMax) {
    let goldengate_time_div = document.getElementById(sum_elem);
    let time_dataset = [];

    for (let i = 0; i < data.length; i++) {
        time_dataset.push(
            {
                label: labels[i],
                data: [data[i]], 
                backgroundColor: time_colors[i]
            }
        );
    }

    const goldengate_time = {
        labels: labels2,
        datasets: time_dataset
    }

    const goldengate_time_config = time_config(goldengate_time, yMax);

    const golengate_time_bar = new Chart(
        document.getElementById(chart_elem), 
        goldengate_time_config
    );

    goldengate_time_div.innerHTML = sum + "hr";
}

function goldengate_cost_chart(data, labels, sum, chart_elem, offset, offset_amt) {
    let [clean_data, clean_labels, clean_colors] = cost_cleanup(data, labels, cost_colors);

    if(offset) {
        clean_colors.push('transparent');
        clean_data.push(offset_amt)
    }
    
    const goldengate_cost = {
        labels: clean_labels,
        datasets: [{
            label: 'My First Dataset',
            data: clean_data,
            backgroundColor: clean_colors,
            hoverOffset: 4,
            borderRadius: 6
        }]
    };

    const goldengate_cost_config = cost_config(goldengate_cost, sum);

    const goldengate_cost_doughnut = new Chart(
        document.getElementById(chart_elem),
        goldengate_cost_config
    );
}

function goldengate_risk_chart(data, labels, chart_elem) {
    const goldengate_risk = {
        labels: labels,
        datasets: [{
            label: 'assembly',
            data: data,
            backgroundColor: risk_colors,
            borderColor: risk_border_colors,
            borderWidth: 2
            // borderWidth: {
            //     top: 4,
            //     bottom: 4
            // }
        }]
    };

    const goldengate_risk_config = risk_config(goldengate_risk);

    const goldengate_risk_bar = new Chart(
        document.getElementById(chart_elem),
        goldengate_risk_config
    );
}