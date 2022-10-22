Chart.defaults.font.size = 15;
Chart.defaults.plugins.title.font.size = 22;
Chart.defaults.plugins.legend.position = 'bottom';
Chart.defaults.plugins.title.align = 'start';

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

function cost_cleanup(data, labels, colors) {
    let new_data = [];
    let new_labels = [];
    let new_colors = [];

    for (let i = 0; i < data.length; i++) {
        if (data[i] > 0) {
            new_labels.push(labels[i]);
            new_data.push(data[i]);
            new_colors.push(colors[i]);
        }
    }

    return [new_data, new_labels, new_colors];
}

function time_config(dataset_obj, yMax) {
    return {
        type: 'bar',
        data: dataset_obj,
        options: {
            maintainAspectRatio: false,
            plugins: {
                title: {
                    display: true,
                    text: 'Time (min)'
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
                    max: yMax,
                    // grid: {
                    //     display: false,    
                    // }
                }
            }
        }
    };
}

function cost_config(dataset_obj, cost_sum) {
    return {
        type: 'doughnut',
        data: dataset_obj,
        options: {
            cutout: '65%',
            plugins: {
                title: {
                    display: true,
                    text: 'Cost ($)'
                },
                counter: {
                    fontColor: '#092210',
                    font_size: 35,
                    font_family: 'sans-serif',
                    cost_sum: cost_sum
                }
            },
        },
        plugins: [counter]
    };
}

function risk_config(dataset_obj) {
    return {
        type: 'bar',
        data: dataset_obj,
        options: {
            maintainAspectRatio: false,
            scales: {
                y: {
                    beginAtZero: true,
                    max: 1,
                    min: -1,
                    title: {
                        display: true,
                        text: 'Log odds of failure'
                    }
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
                },
                subtitle: {
                    display: true,
                    text: '-1: less risk, 1: more risk'
                }
            },
        }
    };
}

// labels
const labels1 = ['Assembly Times'];
const labels2 = ['Assembly Times'];
const labels3 = ['Assembly Times'];

const gibson_time_labels = [];
const goldengate_time_labels = [];
const pcr_time_labels = [];
const slic_time_labels = [];
const biobricks_time_labels = [];

const gibson_cost_labels = ['primers', 'parts', 'genes', 'blocks', 'polymerase'];
const goldengate_cost_labels = ['primers', 'parts', 'genes', 'ligase', 'type2s re'];
const risk_labels = ['thing1', 'thing2', 'thing3', 'thing4'];

// colors 
const time_colors = ['#69130F', '#B72D26', '#FF2A1F'];0
const cost_colors = ['#49EC7A', '#3ACB66', '#2EAA53', '#238A42', '#195C2D','#49EC7A', '#3ACB66', '#2EAA53', '#238A42', '#195C2D', '#49EC7A', '#3ACB66', '#2EAA53', '#238A42', '#195C2D'];
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