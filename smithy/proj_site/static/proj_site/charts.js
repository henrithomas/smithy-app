// // Chartist...
// var data = {
//     labels: ['Bananas - 42%', 'Apples - 25%', 'Grapes - 33%'],
//     series: [5, 3, 4]
// };
// var sum = function(a, b) { return a + b };

// new Chartist.Pie('#barChart', data, {
//     labelInterpolationFnc: function(value) {
//         return value;
//     }
// });

// // Initialize a Line chart in the container with the ID chart1
// new Chartist.Line('#chart1', {
//     labels: [1, 2, 3, 4],
//     series: [[100, 120, 180, 200]]
// });

// // Initialize a Line chart in the container with the ID chart2
// new Chartist.Bar('#chart2', {
//     labels: [1, 2, 3, 4],
//     series: [[5, 2, 8, 3]]
// });


// // Britecharts
// // Instantiate bar chart and container
// const barChart = britecharts.bar();
// const container = d3.select('.bar-container');

// // Create Dataset with proper shape
// const barData = [
//     { name: 'Luminous', value: 2 },
//     { name: 'Glittering', value: 5 },
//     { name: 'Intense', value: 4 },
//     { name: 'Radiant', value: 3 }
// ];

// // Configure chart
// barChart
//     .margin({left: 100})
//     .isHorizontal(true)
//     .height(400)
//     .width(600);

// container.datum(barData).call(barChart); 


let myChart1 = document.getElementById("myChart").getContext("2d");
let labels1 = ['test 1'];
// let data1 = [60, 30, 10];
let colors1 = ['#69130F', '#B72D26', '#FF2A1F'];

let chart1 = new Chart(myChart1, {
    type: 'bar',
    data: {
        labels: labels1,
        datasets: [{
            label: 'pcr',
            data: [30], 
            backgroundColor: '#69130F'
        },
        {
            label: 'digestion',
            data: [15], 
            backgroundColor: '#B72D26'
        },
        {
            label: 'ligation',
            data: [10], 
            backgroundColor: '#FF2A1F'
        }]
    },
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
                max: 100
            }
        }
    }
});

let myChart2 = document.getElementById("myChart2").getContext("2d");
let labels2 = ['test 2'];
// let colors1 = ['#69130F', '#B72D26', '#FF2A1F'];

let chart2 = new Chart(myChart2, {
    type: 'bar',
    data: {
        labels: labels2,
        datasets: [{
            label: 'pcr',
            data: [30], 
            backgroundColor: '#69130F'
        },
        {
            label: 'digestion',
            data: [40], 
            backgroundColor: '#B72D26'
        },
        {
            label: 'ligation',
            data: [10], 
            backgroundColor: '#FF2A1F'
        }]
    },
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
                max: 100
            }
        }
    }
});

let myChart3 = document.getElementById("myChart3").getContext("2d");
let labels3 = ['test 3'];
// let colors1 = ['#69130F', '#B72D26', '#FF2A1F'];

let chart3 = new Chart(myChart3, {
    type: 'bar',
    data: {
        labels: labels3,
        datasets: [{
            label: 'pcr',
            data: [30], 
            backgroundColor: '#69130F'
        },
        {
            label: 'digestion',
            data: [15], 
            backgroundColor: '#B72D26'
        },
        {
            label: 'ligation',
            data: [26], 
            backgroundColor: '#FF2A1F'
        }]
    },
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
                max: 100
            }
        }
    }
});

let myRadar = document.getElementById("myRadar").getContext("2d");
let data2 = [65, 59, 90, 81, 56, 55];
let radLabels = ['primers', 'parts', 'genes', 'blocks', 'type2s re', 'ligase'];

new Chart(myRadar, {
    type: 'doughnut',
    data: {
        labels: radLabels,
        datasets: [{
            label: 'My First Dataset',
            data: [65, 59, 90, 81, 56, 55],
            backgroundColor: [
                'rgb(255, 99, 132)',
                'rgb(54, 162, 235)',
                'rgb(255, 205, 86)'
            ],
            hoverOffset: 4
        }]
    },
    options: {
        plugins: {
            title: {
                display: true,
                text: 'Cost'
            }
        },
    }
});

let myBar = document.getElementById("myBar").getContext("2d");
let data3 = [-.9, .8, 0.1, 0.5, -0.5, .8];
let barLabels = ['primers', 'parts', 'genes', 'blocks', 'type2s re', 'ligase'];

new Chart(myBar, {
    type: 'bar',
    data: {
        labels: barLabels,
        datasets: [{
            label: 'dataset test',
            data: data3,
            backgroundColor: [
            'rgba(255, 99, 132, 0.2)',
            'rgba(255, 159, 64, 0.2)',
            'rgba(255, 205, 86, 0.2)',
            'rgba(75, 192, 192, 0.2)',
            'rgba(54, 162, 235, 0.2)',
            'rgba(153, 102, 255, 0.2)',
            ],
            borderColor: [
            'rgb(255, 99, 132)',
            'rgb(255, 159, 64)',
            'rgb(255, 205, 86)',
            'rgb(75, 192, 192)',
            'rgb(54, 162, 235)',
            'rgb(153, 102, 255)',
            ],
            borderWidth: 1
        }]
      },
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
        },
    }
});