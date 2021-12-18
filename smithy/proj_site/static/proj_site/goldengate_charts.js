function goldengate_time_chart(data, chart_elem, sum_elem) {
    let goldengate_time_div = document.getElementById(sum_elem);
    let time_dataset = [];
    let goldengate_time_sum = arr_sum(data);

    for (let i = 0; i < data.length; i++) {
        time_dataset.push(
            {
                //label: goldengate_time_labels[i]
                label: 'pcr',
                data: [data[i]], 
                backgroundColor: time_colors[i]
            }
        );
    }

    const goldengate_time = {
        labels: labels2,
        datasets: time_dataset
    }

    const goldengate_time_config = {
        type: 'bar',
        data: goldengate_time,
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
                    stacked: true
                },
                y: {
                    stacked: true,
                }
            }
        }
    };

    const golengate_time_bar = new Chart(
        document.getElementById(chart_elem), 
        goldengate_time_config
    );

    goldengate_time_div.innerHTML = goldengate_time_sum + "hr";
}

function goldengate_cost_chart(data, chart_elem) {
    let goldengate_cost_sum = arr_sum(data);

    const goldengate_cost = {
        labels: goldengate_cost_labels,
        datasets: [{
            label: 'My First Dataset',
            data: data,
            backgroundColor: cost_colors,
            hoverOffset: 4,
            borderRadius: 6
        }]
    };

    const goldengate_cost_config = {
        type: 'doughnut',
        data: goldengate_cost,
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
                    cost_sum: goldengate_cost_sum
                }
            },
        },
        plugins: [counter]
    };

    const goldengate_cost_doughnut = new Chart(
        document.getElementById(chart_elem),
        goldengate_cost_config
    );
}

function goldengate_risk_chart(data, chart_elem) {
    const goldengate_risk = {
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

    const goldengate_risk_config = {
        type: 'bar',
        data: goldengate_risk,
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

    const goldengate_risk_bar = new Chart(
        document.getElementById(chart_elem),
        goldengate_risk_config
    );
}