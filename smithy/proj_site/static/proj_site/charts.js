// Chart.js defaults
// Chart.defaults.font.size = 15
// Chart.defaults.plugins.title.font.size = 20
// Chart.defaults.plugins.legend.position = 'bottom'

// checker = !!document.getElementById("gibson-time");
// test_html = document.getElementById("testtext");
// test_html.innerHTML = "TEST";

let gibson_time_div = document.getElementById("gibson-time-sum");
let goldengate_time_div = document.getElementById("goldengate-time-sum");
let pcr_time_div = document.getElementById("pcr-time-sum");

// counter plugin 
const counter = { 
    id: 'counter',
    beforeDraw(chart, args, options) {
        const { ctx, chartArea: { top, right, bottom, left, width, height } } = chart;
        ctx.save();
        
        ctx.font = options.font_size  + 'px ' + options.font_family;
        ctx.textAlign = 'center';
        ctx.fillStyle = options.fontColor;
        ctx.fillText('$' + options.cost_sum, width / 2, top + (height / 2));
    }
} 

function arr_sum(arr) {
    let sum = 0;
    for (let i = 0; i < arr.length; i++) {
        sum += arr[i];
    }
    return sum;
}

function array_max(arr) {
    return Math.max.apply(null, arr);
}

// SETUP
// data
const gibson_times = [27, 15, 12];
const goldengate_times = [30, 42, 9];
const pcr_times = [30, 15, 26];

const gibson_costs = [65, 59, 90, 81, 56, 55];
const goldengate_costs = [65, 59, 90, 81, 161];
const data2 = [65, 59, 90, 81, 56, 55];
let cost_compare = [];

const gibson_risks = [-.9, .8, 0.1, 0.5];

let gibson_time_sum = arr_sum(gibson_times);
let goldengate_time_sum = arr_sum(goldengate_times);
let pcr_time_sum = arr_sum(pcr_times);

let gibson_cost_sum = arr_sum(gibson_costs);
let goldengate_cost_sum = arr_sum(goldengate_costs);

cost_compare.push(gibson_cost_sum);
cost_compare.push(goldengate_cost_sum);
let cost_max = array_max(cost_compare);


// labels
const labels1 = ['Assembly Times'];
const labels2 = ['Assembly Times'];
const labels3 = ['Assembly Times'];
const gibson_cost_labels = ['primers', 'parts', 'genes', 'blocks', 'polymerase'];
const goldengate_cost_labels = ['primers', 'parts', 'genes', 'ligase', 'type2s re'];
const risk_labels = ['thing1', 'thing2', 'thing3', 'thing4'];

// colors 
const time_colors = ['#69130F', '#B72D26', '#FF2A1F'];
// TODO: remove transparent so that it can be appended after the cost differences are calculated
const cost_colors = ['#49EC7A', '#3ACB66', '#2EAA53', '#238A42', '#195C2D', 'transparent'];
// '#0F341A'
const risk_colors = [
    'rgba(92, 230, 230, 0.5)',
    'rgba(74, 185, 185, 0.5)',
    'rgba(59, 146, 146, 0.5)',
    'rgba(51, 125, 125, 0.5)',
    'rgba(14, 41, 41, 0.5)',
    'rgba(29, 73, 73, 0.5)',
    'rgba(19, 49, 49, 0.5)',
    ];
const risk_border_colors = [
    'rgb(92, 230, 230)',
    'rgb(74, 185, 185)',
    'rgb(59, 146, 146)',
    'rgb(51, 125, 125)',
    'rgb(14, 41, 41)',
    'rgb(29, 73, 73)',
    'rgb(19, 49, 49)',
    ]

const time_sum_colors = [];
const cost_sum_colors = [];
const risk_sum_colors = [];

//defaults
Chart.defaults.font.size = 15;
Chart.defaults.plugins.title.font.size = 22;
Chart.defaults.plugins.legend.position = 'bottom';
Chart.defaults.scales.linear.max = 81;
Chart.defaults.plugins.title.align = 'start';

// gibson setup
// const <data_setup> = {
//     labels: <labels>,
//     datasets: [
//         <data>
//     ]
// }

const gibson_time = {
    labels: labels1,
    datasets: [{
        label: 'pcr',
        data: [gibson_times[0]], 
        backgroundColor: time_colors[0]
    },
    {
        label: 'digestion',
        data: [gibson_times[1]], 
        backgroundColor: time_colors[1]
    },
    {
        label: 'ligation',
        data: [gibson_times[2]], 
        backgroundColor: time_colors[2]
    }]
};

const gibson_cost = {
    labels: gibson_cost_labels,
    datasets: [{
        label: 'My First Dataset',
        data: gibson_costs,
        backgroundColor: cost_colors,
        hoverOffset: 4,
        borderRadius: 6
    }]
}

const gibson_risk = {
    labels: risk_labels,
    datasets: [{
        label: 'dataset test',
        data: gibson_risks,
        backgroundColor: risk_colors,
        borderColor: risk_border_colors,
        borderWidth: 2
        // borderWidth: {
        //     top: 4,
        //     bottom: 4
        // }
    }]
}

// goldengate
const goldengate_time = {
    labels: labels2,
    datasets: [{
        label: 'pcr',
        data: [goldengate_times[0]], 
        backgroundColor: time_colors[0]
    },
    {
        label: 'digestion',
        data: [goldengate_times[1]], 
        backgroundColor: time_colors[1]
    },
    {
        label: 'ligation',
        data: [goldengate_times[2]], 
        backgroundColor: time_colors[2]
    }]
};

const goldengate_cost = {
    labels: goldengate_cost_labels,
    datasets: [{
        label: 'My First Dataset',
        data: goldengate_costs,
        backgroundColor: cost_colors,
        hoverOffset: 4,
        borderRadius: 6
    }]
}

// pcr-soe
const pcr_time = {
    labels: labels3,
    datasets: [{
        label: 'pcr',
        data: [pcr_times[0]], 
        backgroundColor: '#69130F'
    },
    {
        label: 'digestion',
        data: [pcr_times[1]], 
        backgroundColor: '#B72D26'
    },
    {
        label: 'ligation',
        data: [pcr_times[2]], 
        backgroundColor: '#FF2A1F'
    }]
}

// slic

// biobricks



// CONFIG
// const config = {
//     type: 'line',
//     data: data,
//     options: {}
//   };

// gibson
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
}

// goldengate
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

// pcr-soe
const pcr_time_config = {
    type: 'bar',
    data: pcr_time,
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
                stacked: true
            }
        }
    }
};

// slic

// biobricks



// RENDER

// gibson
const gibson_time_bar = new Chart(
    document.getElementById("gibson-time"),
    gibson_time_config
);
gibson_time_div.innerHTML = gibson_time_sum + "hr";

const gibson_cost_doughnut = new Chart(
    document.getElementById("gibson-cost"),
    gibson_cost_config
)

const gibson_risk_bar = new Chart(
    document.getElementById("gibson-risk"),
    gibson_risk_config
)

// goldengate
const golengate_time_bar = new Chart(
    document.getElementById("goldengate-time"), 
    goldengate_time_config
);
goldengate_time_div.innerHTML = goldengate_time_sum + "hr";

const goldengate_cost_doughnut = new Chart(
    document.getElementById("goldengate-cost"),
    goldengate_cost_config
);

// pcr-soe
const pcr_time_bar = new Chart(
    document.getElementById("pcr-time"),
    pcr_time_config
);
pcr_time_div.innerHTML = pcr_time_sum + "hr";

// slic

// biobricks
