document.getElementById('poly_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_polymerase_cost').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_polymerase_cost').setAttribute('hidden', true);
    }
});