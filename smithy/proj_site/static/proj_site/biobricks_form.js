document.getElementById('id_EcoRI_cost').addEventListener('change', function() {
    document.getElementById("EcoRI_cost_display").innerHTML = "$" + this.value;
});

document.getElementById('id_XbaI_cost').addEventListener('change', function() {
    document.getElementById("XbaI_cost_display").innerHTML = "$" + this.value;
});

document.getElementById('id_SpeI_cost').addEventListener('change', function() {
    document.getElementById("SpeI_cost_display").innerHTML = "$" + this.value;
});

document.getElementById('id_PstI_cost').addEventListener('change', function() {
    document.getElementById("PstI_cost_display").innerHTML = "$" + this.value;
});