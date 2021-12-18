document.getElementById('id_polymerase_cost').addEventListener('change', function() {
    document.getElementById("poly_cost_display").innerHTML = "$" + this.value;
});