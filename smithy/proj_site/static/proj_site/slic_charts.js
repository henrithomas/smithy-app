function slic_time_chart(data, chart_elem, sum_elem) {
    let slic_time_div = document.getElementById(sum_elem);
    let time_dataset = [];
    let slic_time_sum = arr_sum(data);

    for (let i = 0; i < data.length; i++) {
        time_dataset.push(
            {
                //label: slic_time_labels[i]
                label: 'pcr',
                data: [data[i]], 
                backgroundColor: time_colors[i]
            }
        );
    }

    const slic_time = {
        labels: labels1,
        datasets: time_dataset
    };

    const slic_time_config = time_config(slic_time);

    const slic_time_bar = new Chart(
        document.getElementById(chart_elem),
        slic_time_config
    );

    slic_time_div.innerHTML = slic_time_sum + "hr";
}

function slic_cost_chart(data, chart_elem) {
    let slic_cost_sum = arr_sum(data);

    const slic_cost = {
        labels: slic_cost_labels,
        datasets: [{
            label: 'My First Dataset',
            data: data,
            backgroundColor: cost_colors,
            hoverOffset: 4,
            borderRadius: 6
        }]
    };

    const slic_cost_config = cost_config(slic_cost, slic_cost_sum);

    const slic_cost_doughnut = new Chart(
        document.getElementById(chart_elem),
        slic_cost_config
    );
}

function slic_risk_chart(data, chart_elem) {
    const slic_risk = {
        labels: risk_labels,
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

    const slic_risk_config = risk_config(slic_risk);

    const slic_risk_bar = new Chart(
        document.getElementById(chart_elem),
        slic_risk_config
    );
}