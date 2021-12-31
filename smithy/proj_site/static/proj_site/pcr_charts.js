function pcr_time_chart(data, chart_elem, sum_elem) {
    let pcr_time_div = document.getElementById(sum_elem);
    let time_dataset = [];
    let pcr_time_sum = arr_sum(data);

    for (let i = 0; i < data.length; i++) {
        time_dataset.push(
            {
                //label: gibson_time_labels[i]
                label: 'pcr',
                data: [data[i]], 
                backgroundColor: time_colors[i]
            }
        );
    }

    const pcr_time = {
        labels: labels3,
        datasets: time_dataset
    };

    const pcr_time_config = time_config(pcr_time);

    const pcr_time_bar = new Chart(
        document.getElementById(chart_elem),
        pcr_time_config
    );

    pcr_time_div.innerHTML = pcr_time_sum + "hr";
}

function pcr_cost_chart(data, chart_elem) {
    let pcr_cost_sum = arr_sum(data);

    const pcr_cost = {
        labels: pcr_cost_labels,
        datasets: [{
            label: 'My First Dataset',
            data: data,
            backgroundColor: cost_colors,
            hoverOffset: 4,
            borderRadius: 6
        }]
    };

    const pcr_cost_config = cost_config(pcr_cost, pcr_cost_sum);

    const pcr_cost_doughnut = new Chart(
        document.getElementById(chart_elem),
        pcr_cost_config
    );
}

function pcr_risk_chart(data, chart_elem) {
    const pcr_risk = {
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

    const pcr_risk_config = risk_config(pcr_risk);

    const pcr_risk_bar = new Chart(
        document.getElementById(chart_elem),
        pcr_risk_config
    );
}