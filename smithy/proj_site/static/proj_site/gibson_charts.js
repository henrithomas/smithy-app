function gibson_time_chart(data, chart_elem, sum_elem) {
    let gibson_time_div = document.getElementById(sum_elem);
    let time_dataset = [];
    let gibson_time_sum = arr_sum(data);

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

    const gibson_time = {
        labels: labels1,
        datasets: time_dataset
    };

    const gibson_time_config = {
        type: 'bar',
        data: gibson_time,
        options: {
            maintainAspectRatio: false,
            plugins: {
                title: {
                    display: true,
                    text: 'Time'
                }
            },
            scales: {
                x: {
                    stacked: true,
                    // grid: {
                    //     display: false
                    // }
                },
                y: {
                    stacked: true,
                    // grid: {
                    //     display: false,    
                    // }
                }
            }
        }
    };

    const gibson_time_bar = new Chart(
        document.getElementById(chart_elem),
        gibson_time_config
    );

    gibson_time_div.innerHTML = gibson_time_sum + "hr";
}

function gibson_cost_chart(data, chart_elem) {
    let gibson_cost_sum = arr_sum(data);

    const gibson_cost = {
        labels: gibson_cost_labels,
        datasets: [{
            label: 'My First Dataset',
            data: data,
            backgroundColor: cost_colors,
            hoverOffset: 4,
            borderRadius: 6
        }]
    };

    const gibson_cost_config = {
        type: 'doughnut',
        data: gibson_cost,
        options: {
            cutout: '65%',
            plugins: {
                title: {
                    display: true,
                    text: 'Cost'
                },
                counter: {
                    fontColor: '#092210',
                    font_size: 35,
                    font_family: 'sans-serif',
                    cost_sum: gibson_cost_sum
                }
            },
        },
        plugins: [counter]
    };

    const gibson_cost_doughnut = new Chart(
        document.getElementById(chart_elem),
        gibson_cost_config
    );
}

function gibson_risk_chart(data, chart_elem) {
    const gibson_risk = {
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

    const gibson_risk_config = {
        type: 'bar',
        data: gibson_risk,
        options: {
            maintainAspectRatio: false,
            scales: {
                y: {
                    beginAtZero: true,
                    max: 1,
                    min: -1
                },
            },
            plugins: {
                title: {
                    display: true,
                    text: 'Risk'
                },
                legend: {
                    display: false
                },
                tooltip: {
                    // mode: 'point'
                    enabled: true
                }
            },
        }
    };

    const gibson_risk_bar = new Chart(
        document.getElementById(chart_elem),
        gibson_risk_config
    );
}