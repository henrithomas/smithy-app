function biobricks_time_chart(data, labels, sum, chart_elem, sum_elem) {
    let biobricks_time_div = document.getElementById(sum_elem);
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

    const biobricks_time = {
        labels: labels1,
        datasets: time_dataset
    };

    const biobricks_time_config = time_config(biobricks_time);

    const biobricks_time_bar = new Chart(
        document.getElementById(chart_elem),
        biobricks_time_config
    );

    biobricks_time_div.innerHTML = sum + "hr";
}

function biobricks_cost_chart(data, labels, sum, chart_elem) {

    const biobricks_cost = {
        labels: labels,
        datasets: [{
            label: 'My First Dataset',
            data: data,
            backgroundColor: cost_colors,
            hoverOffset: 4,
            borderRadius: 6
        }]
    };

    const biobricks_cost_config = cost_config(biobricks_cost, sum);

    const biobricks_cost_doughnut = new Chart(
        document.getElementById(chart_elem),
        biobricks_cost_config
    );
}

function biobricks_risk_chart(data, labels, chart_elem) {
    const biobricks_risk = {
        labels: labels,
        datasets: [{
            label: 'dataset test',
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

    const biobricks_risk_config = risk_config(biobricks_risk);

    const biobricks_risk_bar = new Chart(
        document.getElementById(chart_elem),
        biobricks_risk_config
    );
}