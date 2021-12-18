document.getElementById('id_exonuclease_cost').addEventListener('change', function() {
    document.getElementById("exonuc_cost_display").innerHTML = "$" + this.value;
});

document.getElementById('id_ligase_cost').addEventListener('change', function() {
    document.getElementById("ligase_cost_display").innerHTML = "$" + this.value;
});