document.getElementById('id_re_cost').addEventListener('change', function() {
    document.getElementById("re_cost_display").innerHTML = "$" + this.value;
});

document.getElementById('id_ligase_cost').addEventListener('change', function() {
    document.getElementById("ligase_cost_display").innerHTML = "$" + this.value;
});
