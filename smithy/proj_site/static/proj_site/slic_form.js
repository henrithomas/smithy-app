document.getElementById('exo_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_exonuclease_cost').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_exonuclease_cost').setAttribute('hidden', true);
    }
});

document.getElementById('ligase_check').addEventListener('click', function() {
    if (this.checked) {
        document.getElementById('id_ligase_cost').removeAttribute('hidden');
    }
    else {
        document.getElementById('id_ligase_cost').setAttribute('hidden', true);
    }
});