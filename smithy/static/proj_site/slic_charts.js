function slic_time_chart(data, labels, sum, chart_elem, sum_elem, yMax) {
    let slic_time_div = document.getElementById(sum_elem);
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

    const slic_time = {
        labels: labels1,
        datasets: time_dataset
    };

    const slic_time_config = time_config(slic_time, yMax);

    const slic_time_bar = new Chart(
        document.getElementById(chart_elem),
        slic_time_config
    );

    slic_time_div.innerHTML = sum + "hr";
}

function slic_cost_chart(data, labels, sum, chart_elem, offset, offset_amt) {
    let [clean_data, clean_labels, clean_colors] = cost_cleanup(data, labels, cost_colors);

    if(offset) {
        clean_colors.push('transparent');
        clean_data.push(offset_amt)
    }
    
    const slic_cost = {
        labels: clean_labels,
        datasets: [{
            label: 'My First Dataset',
            data: clean_data,
            backgroundColor: clean_colors,
            hoverOffset: 4,
            borderRadius: 6
        }]
    };

    const slic_cost_config = cost_config(slic_cost, sum);

    const slic_cost_doughnut = new Chart(
        document.getElementById(chart_elem),
        slic_cost_config
    );
}

function slic_risk_chart(data, labels, chart_elem) {
    const slic_risk = {
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

    const slic_risk_config = risk_config(slic_risk);

    const slic_risk_bar = new Chart(
        document.getElementById(chart_elem),
        slic_risk_config
    );
}