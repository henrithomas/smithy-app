function gibson_time_chart(data, labels, sum, chart_elem, sum_elem) {
    let gibson_time_div = document.getElementById(sum_elem);
    let time_dataset = [];
    let gibson_time_sum = arr_sum(data);

    for (let i = 0; i < data.length; i++) {
        time_dataset.push(
            {
                //label: gibson_time_labels[i]
                label: labels[i],
                data: [data[i]], 
                backgroundColor: time_colors[i]
            }
        );
    }

    const gibson_time = {
        labels: labels1,
        datasets: time_dataset
    };

    const gibson_time_config = time_config(gibson_time);

    const gibson_time_bar = new Chart(
        document.getElementById(chart_elem),
        gibson_time_config
    );

    gibson_time_div.innerHTML = gibson_time_sum + "hr";
}

function gibson_cost_chart(data, labels, sum, chart_elem) {
    let gibson_cost_sum = arr_sum(data);

    const gibson_cost = {
        labels: labels,
        datasets: [{
            label: 'My First Dataset',
            data: data,
            backgroundColor: cost_colors,
            hoverOffset: 4,
            borderRadius: 6
        }]
    };

    const gibson_cost_config = cost_config(gibson_cost, sum);

    const gibson_cost_doughnut = new Chart(
        document.getElementById(chart_elem),
        gibson_cost_config
    );
}

function gibson_risk_chart(data, labels, chart_elem) {
    const gibson_risk = {
        labels: labels,
        datasets: [{
            label: 'risk',
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

    const gibson_risk_config = risk_config(gibson_risk);

    const gibson_risk_bar = new Chart(
        document.getElementById(chart_elem),
        gibson_risk_config
    );
}